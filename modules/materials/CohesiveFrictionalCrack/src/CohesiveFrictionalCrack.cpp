#include "Marmot/CohesiveFrictionalCrack.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMath.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Eigen;

  CohesiveFrictionalCrack::CohesiveFrictionalCrack( const double* materialProperties,
                                                    int           nMaterialProperties,
                                                    int           materialNumber )
    : Marmot::Materials::MaterialCohesive( materialProperties, nMaterialProperties, materialNumber ),
      // the number of properties is checked before the first one is read
      Kc( nMaterialProperties >= 10
            ? materialProperties[0]
            : throw std::invalid_argument( Marmot::Utilities::MakeString()
                                           << __PRETTY_FUNCTION__ << ": 10 material properties required, "
                                           << nMaterialProperties << " given." ) ),
      ft( materialProperties[1] ),
      GF( materialProperties[2] ),
      Gf( materialProperties[3] ),
      ftInt( materialProperties[4] ),
      beta( materialProperties[5] ),
      kappa( materialProperties[6] ),
      mu( materialProperties[7] ),
      eta( materialProperties[8] ),
      uDilMax( materialProperties[9] ),
      txCrit( nMaterialProperties > 10 ? materialProperties[10] : 0.0 ),
      massPerArea( nMaterialProperties > 11 ? materialProperties[11] : 0.0 ),
      alpha0( ft / Kc ),
      alpha1( 2 * Gf / ft + alpha0 ),
      alphaInt( alpha1 - ( alpha1 - alpha0 ) * ftInt / ft ),
      alphaUlt( 2 * ( GF - Gf ) / ftInt + alphaInt )
  {
    if ( Kc <= 0 || ft <= 0 || Gf <= 0 || beta <= 0 || kappa <= 0 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": K_c, f_t, G_f, beta and kappa must be positive." );
    if ( GF <= Gf )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": G_F must be larger than G_f." );
    if ( ftInt <= 0 || ftInt > ft )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": 0 < f_t,int <= f_t is required." );
    if ( mu < 0 || eta < 0 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": mu and eta must not be negative." );
    if ( eta > 0 && uDilMax <= 0 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": the dilatancy limit must be positive if eta > 0." );
    // a deck of the former layout, with the tangent method at [11] and the mass per area at [12], is rejected
    if ( nMaterialProperties > 12 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__ << ": at most 12 material properties, " << nMaterialProperties
                                   << " given (the mass per area is materialProperties[11])." );
    // a NaN passes every ordering test, so it is rejected explicitly
    if ( !std::isfinite( massPerArea ) || massPerArea < 0 )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": the mass per area must be a finite, non-negative number." );

    initializeStateLayout();
  }

  double CohesiveFrictionalCrack::getMaximumStiffness( const Response&, const Eigen::Vector3d& ) const
  {
    return Kc * ( std::max( 1.0, beta * beta / kappa ) + std::sqrt( mu * eta ) );
  }

  double CohesiveFrictionalCrack::getMassPerArea( const double* ) const
  {
    // as the density of a bulk material, only needed for the inertia
    if ( nMaterialProperties <= 11 )
      throw std::runtime_error( Marmot::Utilities::MakeString()
                                << __PRETTY_FUNCTION__
                                << ": mass per area not specified for this material (materialProperties[11])." );
    return massPerArea;
  }

  void CohesiveFrictionalCrack::computeTraction( Response&        response,
                                                 Tangents&        tangents,
                                                 const Increment& increment ) const
  {
    // the history at the beginning of the increment, read before storeState overwrites it
    const double&         alphaMax = stateLayout.getAs< const double& >( response.stateVars, "alphaMax" );
    const Map< Vector2d > slipOld  = stateLayout.getAs< Map< Vector2d > >( response.stateVars, "slip" );

    // explicit part: traction and tangent for the extrapolated internal variables
    const InternalVariables extrapolated = extrapolateInternalVariables( response.stateVars, increment.dT );
    const TractionResponse  rx           = evaluateTraction( increment.jump, extrapolated, slipOld );

    tangents.dTraction_dJump = rx.dTraction_dJump;

    updateResponse( response, rx, increment.dJump );

    // implicit part: the internal variables for the current jump, the history of the next extrapolation
    storeState( response.stateVars, integrateInternalVariables( increment.jump, alphaMax, slipOld ), increment );
  }

  void CohesiveFrictionalCrack::computeTractionExplicit( Response& response, const Increment& increment ) const
  {
    // the history at the beginning of the increment, read before storeState overwrites it
    const double&         alphaMax = stateLayout.getAs< const double& >( response.stateVars, "alphaMax" );
    const Map< Vector2d > slipOld  = stateLayout.getAs< Map< Vector2d > >( response.stateVars, "slip" );

    // no extrapolation: the traction for the internal variables integrated implicitly for the current jump
    const ImplicitUpdate implicit = integrateInternalVariables( increment.jump, alphaMax, slipOld );

    updateResponse( response,
                    evaluateTraction( increment.jump, implicit.internalVariables, slipOld ),
                    increment.dJump );
    storeState( response.stateVars, implicit, increment );
  }

  CohesiveFrictionalCrack::InternalVariables CohesiveFrictionalCrack::extrapolateInternalVariables(
    const double* stateVars,
    double        dT ) const
  {
    const double& alphaMax          = stateLayout.getAs< const double& >( stateVars, "alphaMax" );
    const double& alphaMaxOld       = stateLayout.getAs< const double& >( stateVars, "alphaMaxOld" );
    const double& slipMultiplier    = stateLayout.getAs< const double& >( stateVars, "slipMultiplier" );
    const double& tangentialJump    = stateLayout.getAs< const double& >( stateVars, "tangentialJump" );
    const double& tangentialJumpOld = stateLayout.getAs< const double& >( stateVars, "tangentialJumpOld" );
    const double& dTOld             = stateLayout.getAs< const double& >( stateVars, "dTOld" );
    const double& contact           = stateLayout.getAs< const double& >( stateVars, "contact" );

    // linear in time, u_T not below 0; nothing is extrapolated in the first increment and for dT = 0
    const double r = dT > 0 && dTOld > 0 ? dT / dTOld : 0.0;
    return { alphaMax + r * ( alphaMax - alphaMaxOld ),
             r * slipMultiplier,
             std::max( tangentialJump + r * ( tangentialJump - tangentialJumpOld ), 0.0 ),
             contact > 0 };
  }

  CohesiveFrictionalCrack::ImplicitUpdate CohesiveFrictionalCrack::integrateInternalVariables(
    const Vector3d& jump,
    double          alphaMaxOld,
    const Vector2d& slipOld ) const
  {
    const double   ux     = jump( 0 );
    const Vector2d uT     = jump.tail< 2 >();
    const double   uTNorm = uT.norm();

    // contact, Eq. 17, with the dilatancy of Eq. 21 and 26
    const double uxContact = ux - computeDilatancy( ux, uTNorm ).uDil;
    const bool   contact   = uxContact <= 0;

    // friction, Eq. 22-23
    const auto [slipMultiplier, slip] = computeSlip( uT, uxContact, slipOld );

    return { { computeAlphaMax( jump, alphaMaxOld ), slipMultiplier, uTNorm, contact }, slip };
  }

  void CohesiveFrictionalCrack::updateResponse( Response&               response,
                                                const TractionResponse& traction,
                                                const Eigen::Vector3d&  dJump )
  {
    // the work of the increment by the trapezoidal rule; what is not stored is dissipated
    const double work = 0.5 * ( response.traction + traction.traction ).dot( dJump );
    response.dissipation += work - ( traction.elasticEnergyDensity - response.elasticEnergyDensity );
    response.elasticEnergyDensity = traction.elasticEnergyDensity;
    response.traction             = traction.traction;
  }

  void CohesiveFrictionalCrack::storeState( double*               stateVars,
                                            const ImplicitUpdate& implicit,
                                            const Increment&      increment ) const
  {
    double&         alphaMax          = stateLayout.getAs< double& >( stateVars, "alphaMax" );
    double&         alphaMaxOld       = stateLayout.getAs< double& >( stateVars, "alphaMaxOld" );
    Map< Vector2d > slip              = stateLayout.getAs< Map< Vector2d > >( stateVars, "slip" );
    double&         slipMultiplier    = stateLayout.getAs< double& >( stateVars, "slipMultiplier" );
    double&         tangentialJump    = stateLayout.getAs< double& >( stateVars, "tangentialJump" );
    double&         tangentialJumpOld = stateLayout.getAs< double& >( stateVars, "tangentialJumpOld" );
    double&         dTOld             = stateLayout.getAs< double& >( stateVars, "dTOld" );
    double&         contact           = stateLayout.getAs< double& >( stateVars, "contact" );
    double&         damage            = stateLayout.getAs< double& >( stateVars, "damage" );
    double&         alphaC            = stateLayout.getAs< double& >( stateVars, "alphaC" );

    const InternalVariables& iv = implicit.internalVariables;

    // a negligible increment, e.g. the remainder that closes a step, changes the history only by round-off, which the
    // next extrapolation would amplify by dT_n+1 / dT: it does not shift the history
    if ( increment.dT > 0 && increment.dT >= minRelativeTimeIncrement * dTOld ) {
      alphaMaxOld       = alphaMax;
      tangentialJumpOld = tangentialJump;
      slipMultiplier    = iv.slipMultiplier;
      dTOld             = increment.dT;
    }
    alphaMax       = iv.alphaMax;
    slip           = implicit.slip;
    tangentialJump = iv.tangentialJump;
    damage         = damageFromAlphaMax( iv.alphaMax );
    contact        = iv.contact ? 1.0 : 0.0;
    alphaC         = std::hypot( std::max( increment.jump( 0 ), 0.0 ), beta / kappa * iv.tangentialJump );
  }

  double CohesiveFrictionalCrack::damageFromAlphaMax( double alphaMax ) const
  {
    if ( alphaMax <= alpha0 )
      return 0.0;
    return std::min( 1.0 - softeningFunction( alphaMax ) / ( Kc * alphaMax ), 1.0 );
  }

  CohesiveFrictionalCrack::TractionResponse CohesiveFrictionalCrack::evaluateTraction(
    const Vector3d&          jump,
    const InternalVariables& internalVariables,
    const Vector2d&          slipOld ) const
  {
    const double&  ux = jump( 0 );
    const Vector2d uT = jump.tail< 2 >();
    const double   bk = beta / kappa;
    const double   d  = damageFromAlphaMax( internalVariables.alphaMax );

    // dilatancy, Eq. 21 and 26, for the given tangential jump; it depends on u_x only through g( t_x )
    const auto dilatancy = computeDilatancy( ux, internalVariables.tangentialJump );

    // contact, Eq. 17
    const double uxContact      = ux - dilatancy.uDil;
    const double dUxContact_dUx = 1. - dilatancy.dUDil_dUx;
    const double H              = internalVariables.contact ? 1.0 : 0.0;

    // friction, Eq. 22-23: the given slip multiplier, applied in the direction of the trial traction
    Vector2d tMu      = Kc * ( uT - slipOld );
    Matrix2d dTMu_dUT = Kc * Matrix2d::Identity();
    // the traction relaxed by the slip of the increment: |t_mu| = |t_tr| - K_c dLambda_mu
    const double slipRelaxation = Kc * internalVariables.slipMultiplier;
    if ( slipRelaxation > 0 ) { // the increment slips (implicit: friction limit reached; Impl-Ex: extrapolated slip)
      const double tTrialNorm = tMu.norm();
      if ( tTrialNorm > slipRelaxation ) {
        const Vector2d m = tMu / tTrialNorm;
        tMu -= slipRelaxation * m;
        dTMu_dUT -= Kc * slipRelaxation / tTrialNorm * ( Matrix2d::Identity() - m * m.transpose() );
      }
      else { // the slip would overtake the jump, which is unphysical: capped at slip = u_T, so t_mu = 0
        tMu.setZero();
        dTMu_dUT.setZero();
      }
    }

    // traction, Eq. 16: t_x = K_c u_x,c ( 1 - d + d H ), t_T = ( 1 - d ) K_c beta^2 / kappa u_T + d t_mu H
    TractionResponse r;
    r.traction( 0 )        = ( 1 - d + d * H ) * Kc * uxContact;
    r.traction.tail< 2 >() = ( 1 - d ) * Kc * beta * bk * uT + d * H * tMu;

    r.dTraction_dJump.setZero();
    r.dTraction_dJump( 0, 0 )                     = ( 1 - d + d * H ) * Kc * dUxContact_dUx;
    r.dTraction_dJump.bottomRightCorner< 2, 2 >() = ( 1 - d ) * Kc * beta * bk * Matrix2d::Identity() +
                                                    d * H * dTMu_dUT;

    // recoverable energy of the springs of Eq. 16
    r.elasticEnergyDensity = 0.5 * ( r.traction( 0 ) * uxContact + ( 1 - d ) * Kc * beta * bk * uT.squaredNorm() +
                                     d * H * tMu.squaredNorm() / Kc );

    return r;
  }

  CohesiveFrictionalCrack::DilatancyResponse CohesiveFrictionalCrack::computeDilatancy( double ux, double uT ) const
  {
    DilatancyResponse r{ 0.0, 0.0 };

    // no dilatancy for eta = 0 or u_T = 0
    if ( eta == 0 || uT <= 0 )
      return r;

    // without the pressure function, or if the faces are not pressed together, g = 1 and du_dil/du_x = 0
    double g      = 1.0;
    double dr_dtx = 1.0;

    if ( txCrit < 0 && ux - uDilMax * std::tanh( eta / uDilMax * uT ) < 0 ) {
      // In contact, t_x is the root of r( t_x ) = t_x - K_c ( u_x - u_dil( g( t_x ), |u_T| ) ), with
      // dr/dt_x = 1 + K_c du_dil/dg dg/dt_x
      const double txRoot = solveContactTraction( ux, uT );
      const auto   p      = pressureFunction( txRoot );

      // The bisection resolves t_x only to round-off, but the traction K_c ( u_x - u_dil( g( t_x ) ) ) differs from the
      // root by the residual, i.e. the error amplified by dr/dt_x >> 1. One Newton step on the branch of g at the root
      // reduces it to second order.
      const double th0     = std::tanh( eta * p.g / uDilMax * uT );
      const double dr0_dtx = 1.0 + Kc * eta * uT * ( 1 - th0 * th0 ) * p.dG_dTx;
      const double tx      = txRoot - ( txRoot - Kc * ( ux - uDilMax * th0 ) ) / dr0_dtx;
      const auto   px      = pressureFunction( tx );

      // the derivative at the corrected root
      const double th = std::tanh( eta * px.g / uDilMax * uT );
      dr_dtx          = 1.0 + Kc * eta * uT * ( 1 - th * th ) * px.dG_dTx;
      g               = px.g;
    }

    r.uDil = uDilMax * std::tanh( eta * g / uDilMax * uT );

    // implicit differentiation: du_dil/du_x = ( dr/dt_x - 1 ) / dr/dt_x
    r.dUDil_dUx = ( dr_dtx - 1.0 ) / dr_dtx;

    return r;
  }

  double CohesiveFrictionalCrack::computeAlphaMax( const Vector3d& jump, double alphaMaxOld ) const
  {
    // equivalent opening, Eq. 18; its history maximum grows only on loading (squares compared, root only if needed)
    const double bk     = beta / kappa;
    const double uxPos  = std::max( jump( 0 ), 0.0 );
    const double alpha2 = uxPos * uxPos + bk * bk * jump.tail< 2 >().squaredNorm();
    return alpha2 > alphaMaxOld * alphaMaxOld ? std::sqrt( alpha2 ) : alphaMaxOld;
  }

  CohesiveFrictionalCrack::SlipResponse CohesiveFrictionalCrack::computeSlip( const Vector2d& uT,
                                                                              double          uxContact,
                                                                              const Vector2d& slipOld ) const
  {
    // Coulomb friction limit mu |t_x| with the contact pressure t_x = K_c u_x,c; out of contact it vanishes
    const double frictionLimit = mu * Kc * std::max( -uxContact, 0.0 );

    // elastic (stick) predictor and radial return onto the friction limit
    const Vector2d tTrial      = Kc * ( uT - slipOld );
    const double   tTrialNorm2 = tTrial.squaredNorm();
    if ( tTrialNorm2 <= frictionLimit * frictionLimit )
      return { 0.0, slipOld };

    // slip ( dLambda > 0 ) requires |t_mu| = |t_tr| - K_c dLambda = mu |t_x| (Kuhn-Tucker), which gives dLambda in
    // closed form
    const double tTrialNorm     = std::sqrt( tTrialNorm2 );
    const double slipMultiplier = ( tTrialNorm - frictionLimit ) / Kc;
    return { slipMultiplier, slipOld + slipMultiplier * ( tTrial / tTrialNorm ) };
  }

  double CohesiveFrictionalCrack::solveContactTraction( double ux, double uT ) const
  {
    // The residual r( t_x ) = t_x - K_c ( u_x - u_dil( g( t_x ) ) ) increases monotonically, so the root is unique.
    // (Taking g from the traction of the previous increment instead is unstable: its gain is K_c eta u_T / |t_x,crit|.)
    if ( !std::isfinite( ux ) || !std::isfinite( uT ) )
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << "CohesiveFrictionalCrack: non-finite jump, u_x = " << ux << ", |u_T| = " << uT );

    const auto residual = [&]( double tx ) {
      return tx - Kc * ( ux - uDilMax * std::tanh( eta * pressureFunction( tx ).g / uDilMax * uT ) );
    };

    const double uDil1 = uDilMax * std::tanh( eta / uDilMax * uT );
    // bracket: lo is the traction with g = 1 (full dilatancy), hi the one with g = 0, but at most 0 (contact)
    double lo = Kc * ( ux - uDil1 );
    double hi = std::min( 0.0, Kc * ux );

    // t_x = K_c ( u_x - u_dil ) is resolvable only to the round-off of this difference; on that scale the bracket
    // needs at most log2( 1e14 ) < 47 halvings
    const double tolerance = 1e-14 * Kc * ( std::abs( ux ) + uDil1 );
    for ( int i = 0; hi - lo > tolerance; i++ ) {
      if ( i == 100 )
        throw std::runtime_error( Marmot::Utilities::MakeString()
                                  << "CohesiveFrictionalCrack: the contact traction did not converge, u_x = " << ux
                                  << ", |u_T| = " << uT << ", bracket [" << lo << ", " << hi << "]" );
      const double mid                  = 0.5 * ( lo + hi );
      ( residual( mid ) > 0 ? hi : lo ) = mid;
    }
    return 0.5 * ( lo + hi );
  }

  CohesiveFrictionalCrack::PressureFunctionResponse CohesiveFrictionalCrack::pressureFunction( double tx ) const
  {
    if ( txCrit >= 0 || tx >= 0 )
      return { 1.0, 0.0 };
    if ( tx <= txCrit )
      return { 0.0, 0.0 };
    return { ( txCrit - tx ) / txCrit, -1.0 / txCrit };
  }

  double CohesiveFrictionalCrack::softeningFunction( double alpha ) const
  {
    if ( alpha <= alpha0 )
      return Kc * alpha;
    if ( alpha <= alphaInt )
      return ft * ( alpha1 - alpha ) / ( alpha1 - alpha0 );
    if ( alpha <= alphaUlt )
      return ftInt * ( alphaUlt - alpha ) / ( alphaUlt - alphaInt );
    return 0.0;
  }

} // namespace Marmot::Materials
