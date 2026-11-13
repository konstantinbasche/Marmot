#include "Marmot/CohesiveFrictionalCrack.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialCohesiveFactory.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace Marmot::Testing;
using Marmot::Materials::CohesiveFrictionalCrack;

namespace {

  // concrete properties of Gudzulic & Meschke (2025), Table 2, in N and mm
  std::vector< double > paperProperties()
  {
    //       Kc     ft   GF    Gf     ftInt beta kappa mu   eta  uDilMax
    return { 3.4e7, 2.2, 0.14, 0.056, 0.55, 2.0, 10.0, 0.7, 0.2, 16.0 };
  }

  /// Update of a material point: Impl-Ex for implicit solvers, or the explicit update.
  enum class Update { implEx, explicitDynamics };

  /// A material point that keeps the history of a cohesive law between increments.
  struct CohesivePoint {
    std::vector< double >                                  properties;
    Update                                                 updateType;
    std::unique_ptr< Marmot::Materials::MaterialCohesive > material;
    std::vector< double >                                  stateVars;
    Eigen::Vector3d                                        traction             = Eigen::Vector3d::Zero();
    Eigen::Vector3d                                        jump                 = Eigen::Vector3d::Zero();
    Eigen::Matrix3d                                        tangent              = Eigen::Matrix3d::Zero();
    double                                                 elasticEnergyDensity = 0;
    double                                                 dissipation          = 0;

    double dT = 1.0; ///< time increment of the next update

    explicit CohesivePoint( std::vector< double > props, Update update = Update::implEx )
      : properties( std::move( props ) ), updateType( update )
    {
      material = std::unique_ptr< Marmot::Materials::MaterialCohesive >(
        Marmot::Registry::MaterialCohesiveFactory::createMaterial( "COHESIVEFRICTIONALCRACK",
                                                                   properties.data(),
                                                                   properties.size(),
                                                                   1 ) );
      stateVars.resize( material->getNumberOfRequiredStateVars() );
      material->initializeYourself( stateVars.data(), stateVars.size() );
    }

    /// Update for a trial jump from the current history, which is kept unless @p accept.
    Marmot::Materials::MaterialCohesive::Response update( const Eigen::Vector3d&                         newJump,
                                                          Marmot::Materials::MaterialCohesive::Tangents& tangents,
                                                          bool                                           accept )
    {
      std::vector< double >                                stateVarsCopy = stateVars;
      Marmot::Materials::MaterialCohesive::Response        response{ traction,
                                                              accept ? stateVars.data() : stateVarsCopy.data(),
                                                              elasticEnergyDensity,
                                                              dissipation };
      const Marmot::Materials::MaterialCohesive::Increment increment{ newJump, newJump - jump, 0.0, dT };
      if ( updateType == Update::implEx )
        material->computeTraction( response, tangents, increment );
      else
        material->computeTractionExplicit( response, increment );
      return response;
    }

    /// Advance to a new total jump.
    void moveTo( const Eigen::Vector3d& newJump )
    {
      Marmot::Materials::MaterialCohesive::Tangents tangents;
      const auto                                    response = update( newJump, tangents, true );
      traction                                               = response.traction;
      elasticEnergyDensity                                   = response.elasticEnergyDensity;
      dissipation                                            = response.dissipation;
      tangent                                                = tangents.dTraction_dJump;
      jump                                                   = newJump;
    }

    /// Traction for a trial jump, without changing the history.
    Eigen::Vector3d trialTraction( const Eigen::Vector3d& trialJump )
    {
      Marmot::Materials::MaterialCohesive::Tangents tangents;
      return update( trialJump, tangents, false ).traction;
    }

    /// A point with the same history, updated by @p update.
    CohesivePoint copy( Update update ) const
    {
      CohesivePoint c( properties, update );
      c.stateVars            = stateVars;
      c.traction             = traction;
      c.jump                 = jump;
      c.elasticEnergyDensity = elasticEnergyDensity;
      c.dissipation          = dissipation;
      c.dT                   = dT;
      return c;
    }

    double getState( const std::string& name, int i = 0 )
    {
      return material->getStateView( name, stateVars.data() ).stateLocation[i];
    }

    const CohesiveFrictionalCrack& model() const { return dynamic_cast< const CohesiveFrictionalCrack& >( *material ); }
  };

  /// Load monotonically along a straight path from the current jump to the target.
  void loadTo( CohesivePoint& point, const Eigen::Vector3d& target, int nSteps )
  {
    const Eigen::Vector3d start = point.jump;
    for ( int i = 1; i <= nSteps; i++ )
      point.moveTo( start + ( target - start ) * double( i ) / nSteps );
  }

  /// Compare the Impl-Ex tangent of the next increment with central differences.
  void checkTangent( CohesivePoint& point, const Eigen::Vector3d& nextJump, const std::string& where )
  {
    CohesivePoint probe = point.copy( Update::implEx );

    // the tangent of the increment towards nextJump, from the same history
    Marmot::Materials::MaterialCohesive::Tangents tangents;
    point.update( nextJump, tangents, false );
    const Eigen::Matrix3d C = tangents.dTraction_dJump;

    const double    h = 1e-6 * nextJump.norm() + 1e-14;
    Eigen::Matrix3d CNum;
    for ( int j = 0; j < 3; j++ ) {
      Eigen::Vector3d dJ = Eigen::Vector3d::Zero();
      dJ( j )            = h;
      CNum.col( j )      = ( probe.trialTraction( nextJump + dJ ) - probe.trialTraction( nextJump - dJ ) ) / ( 2 * h );
    }

    const double scale = point.model().Kc;
    throwExceptionOnFailure( ( C - CNum ).norm() <= 1e-5 * scale,
                             Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__ << " (" << where << "): tangent\n"
                                                             << C << "\ndiffers from numerical tangent\n"
                                                             << CNum );

    // with the history fixed by the extrapolation, the tangent is symmetric and positive semi-definite
    throwExceptionOnFailure( ( C - C.transpose() ).norm() <= 1e-12 * scale,
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " (" << where << "): unsymmetric tangent\n"
                               << C );
    const double minEigenvalue = Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d >( C ).eigenvalues().minCoeff();
    throwExceptionOnFailure( minEigenvalue >= -1e-12 * scale,
                             Marmot::Utilities::MakeString()
                               << __PRETTY_FUNCTION__ << " (" << where << "): tangent not positive "
                               << "semi-definite, eigenvalue " << minEigenvalue << "\n"
                               << C );
  }

  /// The damage for a given history maximum, Eq. 19-20.
  double damageOf( const CohesiveFrictionalCrack& par, double alphaMax )
  {
    if ( alphaMax <= par.alpha0 )
      return 0.0;
    const double q = alphaMax <= par.alphaInt ? par.ft * ( par.alpha1 - alphaMax ) / ( par.alpha1 - par.alpha0 )
                     : alphaMax <= par.alphaUlt
                       ? par.ftInt * ( par.alphaUlt - alphaMax ) / ( par.alphaUlt - par.alphaInt )
                       : 0.0;
    return std::min( 1.0 - q / ( par.Kc * alphaMax ), 1.0 );
  }

  /// Exposes the submodels of the material for testing.
  struct CohesiveFrictionalCrackSubmodels : CohesiveFrictionalCrack {
    using CohesiveFrictionalCrack::CohesiveFrictionalCrack;
    using CohesiveFrictionalCrack::computeDilatancy;
    using CohesiveFrictionalCrack::pressureFunction;
    using CohesiveFrictionalCrack::solveContactTraction;
  };

} // namespace

/* ------------------------------------------------------------------------------------------------------------------ *
 * The law, by the explicit update: the traction of the implicitly integrated internal variables for the current jump *
 * ------------------------------------------------------------------------------------------------------------------ */

void testExplicitModeIPeakAndEnergy()
{
  CohesivePoint p( paperProperties(), Update::explicitDynamics );
  const auto&   par = p.model();

  const int    nSteps = 200000;
  const double uxEnd  = 1.2 * par.alphaUlt;

  double peak = 0, energy = 0;
  for ( int i = 1; i <= nSteps; i++ ) {
    const double tOld = p.traction( 0 );
    const double uOld = p.jump( 0 );
    p.moveTo( Eigen::Vector3d( uxEnd * i / nSteps, 0, 0 ) );
    peak = std::max( peak, p.traction( 0 ) );
    energy += 0.5 * ( tOld + p.traction( 0 ) ) * ( p.jump( 0 ) - uOld );
  }

  // the area under the bilinear law of Eq. 19, including the (recoverable) elastic branch
  const double area = 0.5 * par.ft * par.alpha0 + par.GF - 0.5 * par.ftInt * ( par.alpha1 - par.alphaInt );

  throwExceptionOnFailure( checkIfEqual( peak, par.ft, 1e-4 * par.ft ),
                           Marmot::Utilities::MakeString() << "peak " << peak );
  throwExceptionOnFailure( checkIfEqual( energy, area, 1e-4 * area ),
                           Marmot::Utilities::MakeString() << "energy " << energy << " vs " << area );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.0, 1e-12 ), "traction after full separation" );
  throwExceptionOnFailure( checkIfEqual( p.getState( "damage" ), 1.0, 1e-14 ), "damage after full separation" );

  // after complete separation, all the work has been dissipated
  throwExceptionOnFailure( checkIfEqual( p.elasticEnergyDensity, 0.0, 1e-14 ), "stored energy after full separation" );
  throwExceptionOnFailure( checkIfEqual( p.dissipation, energy, 1e-12 * energy ),
                           Marmot::Utilities::MakeString()
                             << "dissipation " << p.dissipation << " vs work " << energy );
}

void testExplicitEnergies()
{
  // without dilatancy the law is a damaged spring with friction: the elastic range stores the work exactly
  auto noDilatancy = paperProperties();
  noDilatancy[8]   = 0.0;

  CohesivePoint p( noDilatancy, Update::explicitDynamics );
  const auto&   par = p.model();

  const Eigen::Vector3d elastic = Eigen::Vector3d( 0.6, 0.3, -0.2 ) * par.alpha0;
  p.moveTo( elastic );
  const double psiElastic = 0.5 * par.Kc *
                            ( std::pow( elastic( 0 ), 2 ) +
                              par.beta * par.beta / par.kappa * elastic.tail< 2 >().squaredNorm() );
  throwExceptionOnFailure( checkIfEqual( p.elasticEnergyDensity, psiElastic, 1e-12 * psiElastic ), "elastic energy" );
  throwExceptionOnFailure( std::abs( p.dissipation ) <= 1e-12 * psiElastic, "no dissipation in the elastic range" );

  // softening dissipates; secant unloading returns the stored energy and dissipates nothing
  loadTo( p, Eigen::Vector3d( 3, 1, 0 ) * par.alpha0 * 10, 50 );
  const double dissipationSoftened = p.dissipation;
  throwExceptionOnFailure( dissipationSoftened > 0, "softening must dissipate" );
  loadTo( p, Eigen::Vector3d( 3, 1, 0 ) * par.alpha0 * 2, 5 );
  throwExceptionOnFailure( checkIfEqual( p.dissipation, dissipationSoftened, 1e-12 * dissipationSoftened ),
                           "unloading must not dissipate" );

  // The dilatancy is not derived from a potential, already in the elastic range: t_x = K_c ( u_x - u_dil( u_T ) )
  // depends on u_T, but t_T not on u_x. A closed elastic cycle ( u_x, 0 ) -> ( u_x, u_T ) -> ( 0, u_T ) -> 0 does the
  // net work K_c u_dil u_x, which ends up as dissipation, with nothing stored.
  {
    CohesivePoint c( paperProperties(), Update::explicitDynamics );
    const auto&   cp = c.model();
    const double  ux = 0.5 * cp.alpha0, uT = 0.5 * cp.alpha0;
    const double  uDil = cp.uDilMax * std::tanh( cp.eta / cp.uDilMax * uT );
    for ( const auto& corner : { Eigen::Vector3d( ux, 0, 0 ),
                                 Eigen::Vector3d( ux, uT, 0 ),
                                 Eigen::Vector3d( 0, uT, 0 ),
                                 Eigen::Vector3d( 0, 0, 0 ) } )
      loadTo( c, corner, 4 );
    throwExceptionOnFailure( checkIfEqual( c.getState( "damage" ), 0.0, 0.0 ), "the cycle must stay elastic" );
    throwExceptionOnFailure( checkIfEqual( c.elasticEnergyDensity, 0.0, 1e-20 ), "nothing stored after the cycle" );
    throwExceptionOnFailure( checkIfEqual( c.dissipation, cp.Kc * uDil * ux, 1e-9 * cp.Kc * uDil * ux ),
                             Marmot::Utilities::MakeString()
                               << "work of the elastic cycle " << c.dissipation << " vs " << cp.Kc * uDil * ux );
  }

  // sliding of a separated crack at constant pressure dissipates the frictional work
  CohesivePoint s( noDilatancy, Update::explicitDynamics );
  loadTo( s, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  const double uxClosed = -1e-7;
  loadTo( s, Eigen::Vector3d( uxClosed, 0, 0 ), 10 );
  loadTo( s, Eigen::Vector3d( uxClosed, 0.01, 0 ), 10 );
  const double dissipationSliding = s.dissipation;
  loadTo( s, Eigen::Vector3d( uxClosed, 0.03, 0 ), 10 );
  const double frictionWork = par.mu * par.Kc * -uxClosed * 0.02;
  throwExceptionOnFailure( checkIfEqual( s.dissipation - dissipationSliding, frictionWork, 1e-10 * frictionWork ),
                           Marmot::Utilities::MakeString()
                             << "frictional dissipation " << s.dissipation - dissipationSliding << " vs "
                             << frictionWork );

  // The dilatancy is not derived from a potential. Shearing a separated crack at a fixed normal jump of zero, the
  // faces slide at once against the pressure K_c eta |u_T| and D grows as K_c eta / 2 ( mu - eta - mu^2 eta ) |u_T|^2
  // while tanh is linear: it decreases for eta > mu / ( 1 + mu^2 ).
  for ( const double mu : { 0.7, 0.1 } ) {
    auto props = paperProperties();
    props[7]   = mu;
    CohesivePoint d( props, Update::explicitDynamics );
    loadTo( d, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
    loadTo( d, Eigen::Vector3d( 0, 0, 0 ), 10 );
    const double dissipationClosed = d.dissipation;
    const double uT                = 1e-3;
    loadTo( d, Eigen::Vector3d( 0, uT, 0 ), 10 );
    const auto&  dp       = d.model();
    const double expected = 0.5 * dp.Kc * dp.eta * ( mu - dp.eta - mu * mu * dp.eta ) * uT * uT;
    throwExceptionOnFailure( checkIfEqual( d.dissipation - dissipationClosed, expected, 1e-6 * std::abs( expected ) ),
                             Marmot::Utilities::MakeString()
                               << "dissipation of dilatant sliding with mu = " << mu << ": "
                               << d.dissipation - dissipationClosed << " vs " << expected );
  }
}

void testExplicitSofteningBranchesAndUnloading()
{
  CohesivePoint p( paperProperties(), Update::explicitDynamics );
  const auto&   par = p.model();

  // on the first branch, the traction follows Eq. 19
  const double ux1 = 0.5 * ( par.alpha0 + par.alphaInt );
  loadTo( p, Eigen::Vector3d( ux1, 0, 0 ), 100 );
  const double q1 = par.ft * ( par.alpha1 - ux1 ) / ( par.alpha1 - par.alpha0 );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), q1, 1e-10 ), "first softening branch" );

  // unloading follows the secant to the origin at once, with frozen damage
  const double d = p.getState( "damage" );
  p.moveTo( Eigen::Vector3d( 0.5 * ux1, 0, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.5 * q1, 1e-10 ), "secant unloading" );
  throwExceptionOnFailure( checkIfEqual( p.getState( "damage" ), d, 1e-15 ), "damage frozen during unloading" );

  // on the second branch
  const double ux2 = 0.5 * ( par.alphaInt + par.alphaUlt );
  loadTo( p, Eigen::Vector3d( ux2, 0, 0 ), 100 );
  const double q2 = par.ftInt * ( par.alphaUlt - ux2 ) / ( par.alphaUlt - par.alphaInt );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), q2, 1e-10 ), "second softening branch" );
}

void testExplicitFrictionalSliding()
{
  auto props = paperProperties();
  props[8]   = 0.0; // no dilatancy

  CohesivePoint p( props, Update::explicitDynamics );
  const auto&   par = p.model();

  // open the crack completely, then close it under compression and shear it
  loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  const double uxClosed = -1e-7;
  loadTo( p, Eigen::Vector3d( uxClosed, 0, 0 ), 100 );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), par.Kc * uxClosed, 1e-9 ), "penalty contact" );

  loadTo( p, Eigen::Vector3d( uxClosed, 0.3, 0.4 ), 100 );
  const double tT = p.traction.tail< 2 >().norm();
  throwExceptionOnFailure( checkIfEqual( tT, par.mu * par.Kc * -uxClosed, 1e-9 ),
                           Marmot::Utilities::MakeString() << "sliding " << tT );

  // the frictional traction points in the sliding direction
  throwExceptionOnFailure( checkIfEqual( p.traction( 2 ) / p.traction( 1 ), 0.4 / 0.3, 1e-10 ), "sliding direction" );
  throwExceptionOnFailure( checkIfEqual( p.getState( "slip", 1 ) / p.getState( "slip", 0 ), 0.4 / 0.3, 1e-10 ),
                           "slip direction" );

  // reversing the shear first sticks, then slides in the opposite direction
  p.moveTo( Eigen::Vector3d( uxClosed, 0.3 - 1e-9, 0.4 ) );
  throwExceptionOnFailure( p.traction( 1 ) < par.mu * par.Kc * -uxClosed * 0.6, "elastic stick on reversal" );
  loadTo( p, Eigen::Vector3d( uxClosed, -0.3, -0.4 ), 100 );
  throwExceptionOnFailure( checkIfEqual( p.traction.tail< 2 >().norm(), par.mu * par.Kc * -uxClosed, 1e-9 ) &&
                             p.traction( 1 ) < 0,
                           "reverse sliding" );
}

void testExplicitDilatancy()
{
  CohesivePoint p( paperProperties(), Update::explicitDynamics );
  const auto&   par = p.model();

  loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  loadTo( p, Eigen::Vector3d( 0, 0, 0 ), 100 );

  // shearing a cracked, normally restrained interface produces compression from the prevented dilatancy, at once
  const double uy = 2.0;
  p.moveTo( Eigen::Vector3d( 0, uy, 0 ) );
  const double uDil = par.uDilMax * std::tanh( par.eta / par.uDilMax * uy );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), -par.Kc * uDil, 1e-9 * par.Kc * uDil ),
                           "dilatancy pressure" );

  // with free normal opening along the dilatancy curve, the faces are just not in contact
  p.moveTo( Eigen::Vector3d( uDil, uy, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction.norm(), 0.0, 1e-9 ), "free dilatant opening" );
}

void testExplicitDilatancyPressureCutoff()
{
  auto props = paperProperties();
  props.push_back( -1.0 ); // t_x,crit

  CohesivePoint p( props, Update::explicitDynamics );
  const auto&   par = p.model();

  loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  loadTo( p, Eigen::Vector3d( -1e-6, 0, 0 ), 100 ); // normal traction -34 < t_x,crit
  p.moveTo( Eigen::Vector3d( -1e-6, 0.01, 0 ) );

  // t_x = -34 < t_x,crit: g = 0, no dilatancy
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), -par.Kc * 1e-6, 1e-9 ), "dilatancy pressure cutoff" );

  // between t_x,crit and 0, g( t_x ) is evaluated with the normal traction it produces: for u_x = 0 and a small
  // tangential jump, t_x = -a g = -a ( 1 + t_x ) with a = K_c u_dil(g = 1)
  CohesivePoint fresh( props, Update::explicitDynamics );
  const double  uT = 0.5 / ( par.Kc * par.eta );
  fresh.moveTo( Eigen::Vector3d( 0, uT, 0 ) );
  const double a = par.Kc * par.uDilMax * std::tanh( par.eta / par.uDilMax * uT );
  throwExceptionOnFailure( checkIfEqual( fresh.traction( 0 ), -a / ( 1 + a ), 1e-9 ),
                           Marmot::Utilities::MakeString()
                             << "implicit pressure function: " << fresh.traction( 0 ) << " vs " << -a / ( 1 + a ) );
}

/// The explicit update neither extrapolates nor depends on the time increment.
void testExplicitDoesNotExtrapolate()
{
  CohesivePoint p( paperProperties() );
  const auto&   par = p.model();
  p.moveTo( Eigen::Vector3d( 0.5 * par.alpha0, 0, 0 ) );

  // a jump beyond the strength: Impl-Ex extrapolates the elastic history, the explicit update softens at once
  const Eigen::Vector3d next( 3 * par.alpha0, 0, 0 );
  const double          softened = par.ft * ( par.alpha1 - next( 0 ) ) / ( par.alpha1 - par.alpha0 );
  for ( const double dT : { 1.0, 5.0, 1e-3 } ) {
    auto e = p.copy( Update::explicitDynamics );
    e.dT   = dT;
    e.moveTo( next );
    throwExceptionOnFailure( checkIfEqual( e.traction( 0 ), softened, 1e-10 * par.ft ),
                             Marmot::Utilities::MakeString() << "explicit traction with dT = " << dT << ": "
                                                             << e.traction( 0 ) << " vs " << softened );
  }

  auto i = p.copy( Update::implEx );
  i.moveTo( next );
  throwExceptionOnFailure( checkIfEqual( i.traction( 0 ), par.Kc * next( 0 ), 1e-10 * par.ft ),
                           Marmot::Utilities::MakeString() << "Impl-Ex traction " << i.traction( 0 ) );
}

/// Both updates store the same implicit integration of the internal variables: along any path, the state variables
/// are identical.
void testImplExAndExplicitStates()
{
  CohesivePoint implex( paperProperties() );
  CohesivePoint explicitPoint( paperProperties(), Update::explicitDynamics );
  const auto&   par = implex.model();

  const std::vector< std::tuple< Eigen::Vector3d, int > >
    path{ { Eigen::Vector3d( 2.0, 1.0, 0.5 ) * par.alpha0, 20 },        // mixed-mode softening
          { Eigen::Vector3d( 0.3 * par.alphaUlt, 0.01, 0.0 ), 30 },     // second branch
          { Eigen::Vector3d( 0.1 * par.alphaUlt, 0.005, 0.0 ), 10 },    // unloading
          { Eigen::Vector3d( -1e-6, 0.0, 0.0 ), 10 },                   // closing
          { Eigen::Vector3d( -1e-6, 0.05, 0.02 ), 20 },                 // sliding with dilatancy
          { Eigen::Vector3d( 1.2 * par.alphaUlt, -0.02, 0.01 ), 40 } }; // separation

  for ( const auto& [target, nSteps] : path ) {
    const Eigen::Vector3d start = implex.jump;
    for ( int i = 1; i <= nSteps; i++ ) {
      const Eigen::Vector3d jump = start + ( target - start ) * double( i ) / nSteps;
      implex.moveTo( jump );
      explicitPoint.moveTo( jump );
      throwExceptionOnFailure( implex.stateVars == explicitPoint.stateVars,
                               Marmot::Utilities::MakeString()
                                 << "the states of the Impl-Ex and the explicit update differ at jump "
                                 << jump.transpose() );
    }
  }
}

void testShearStrength()
{
  // without dilatancy and friction the failure surface is the one of Camacho & Ortiz, Eq. 25:
  // the shear strength is beta * f_t
  auto props = paperProperties();
  props[7]   = 0.0;
  props[8]   = 0.0;

  for ( const auto update : { Update::implEx, Update::explicitDynamics } ) {
    CohesivePoint p( props, update );
    const auto&   par = p.model();

    const int    nSteps = 20000;
    const double uyEnd  = 3 * par.kappa / par.beta * par.alpha0;
    double       peak   = 0;
    for ( int i = 1; i <= nSteps; i++ ) {
      p.moveTo( Eigen::Vector3d( 0, uyEnd * i / nSteps, 0 ) );
      peak = std::max( peak, p.traction( 1 ) );
    }

    throwExceptionOnFailure( checkIfEqual( peak, par.beta * par.ft, 1e-3 * par.ft ),
                             Marmot::Utilities::MakeString() << "shear peak " << peak );
    throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.0, 1e-12 ), "no normal traction in pure shear" );
  }
}

void testEquivalentOpening()
{
  for ( const auto update : { Update::implEx, Update::explicitDynamics } ) {
    CohesivePoint p( paperProperties(), update );
    const auto&   par = p.model();
    const double  bk  = par.beta / par.kappa;

    // on loading, the current equivalent opening is the history maximum
    loadTo( p, Eigen::Vector3d( 0.01, 0.003, -0.004 ), 50 );
    const double alpha = std::hypot( 0.01, bk * 0.005 );
    throwExceptionOnFailure( checkIfEqual( p.getState( "alphaC" ), alpha, 1e-12 ), "alpha_c on loading" );
    throwExceptionOnFailure( checkIfEqual( p.getState( "alphaMax" ), alpha, 1e-12 ), "alpha_max on loading" );

    // closed: only the tangential jump contributes, and the history maximum is kept
    loadTo( p, Eigen::Vector3d( -1e-4, 0.001, 0 ), 20 );
    throwExceptionOnFailure( checkIfEqual( p.getState( "alphaC" ), bk * 0.001, 1e-12 ), "alpha_c closed" );
    throwExceptionOnFailure( checkIfEqual( p.getState( "alphaMax" ), alpha, 1e-12 ), "alpha_max kept" );
  }
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * The Impl-Ex update for implicit solvers                                                                            *
 * ------------------------------------------------------------------------------------------------------------------ */

void testSubmodelDerivatives()
{
  // the derivative of the dilatancy with respect to the normal jump, which enters the Impl-Ex tangent, with g = 1, and
  // with the pressure function: implicit for t_x,crit < t_x < 0, and g = 0 below
  const CohesiveFrictionalCrackSubmodels m( paperProperties().data(), 10, 1 );
  auto                                   cutoffProperties = paperProperties();
  cutoffProperties[9]                                     = 0.1;
  cutoffProperties.push_back( -3.0 );
  const CohesiveFrictionalCrackSubmodels mc( cutoffProperties.data(), 11, 1 );

  const Eigen::Vector2d uT( 0.04, 0.03 );
  const double          uDilHalf = mc.uDilMax * std::tanh( mc.eta * 0.5 / mc.uDilMax * uT.norm() );
  for ( const auto& [model, ux, h, what] :
        std::vector< std::tuple< const CohesiveFrictionalCrackSubmodels*,
                                 double,
                                 double,
                                 std::string > >{ { &m, 0.0, 1e-7, "dilatancy, g = 1" },
                                                  { &mc, 0.02, 1e-7, "dilatancy with pressure function, open" },
                                                  { &mc,
                                                    uDilHalf - 1.5 / mc.Kc,
                                                    1e-12,
                                                    "dilatancy with pressure function, t_x = -1.5 (g = 0.5)" },
                                                  { &mc,
                                                    -1e-3,
                                                    1e-7,
                                                    "dilatancy with pressure function, t_x < t_x,crit (g = 0)" } } ) {
    const double analytical = model->computeDilatancy( ux, uT.norm() ).dUDil_dUx;
    const double numerical  = ( model->computeDilatancy( ux + h, uT.norm() ).uDil -
                               model->computeDilatancy( ux - h, uT.norm() ).uDil ) /
                             ( 2 * h );
    throwExceptionOnFailure( std::abs( analytical - numerical ) <= 1e-5 * std::max( std::abs( numerical ), 1.0 ),
                             Marmot::Utilities::MakeString() << what << ": derivative " << analytical
                                                             << " differs from central differences " << numerical );
  }
}

void testContactTraction()
{
  // the normal traction in contact with the pressure function, i.e. the root of
  // r( t_x ) = t_x - K_c ( u_x - u_dil( g( t_x ) ) ), by bisection: K_c u_x on the plateau g = 0, an error at
  // round-off, monotonic in u_x, and an exception for a non-finite jump
  for ( const auto& [kc, uDilMax, txCrit] :
        std::vector< std::tuple< double, double, double > >{ { 3.4e7, 16.0, -3.0 },
                                                             { 3.4e5, 0.05, -0.01 },
                                                             { 3.4e7, 0.05, -100.0 } } ) {
    auto props = paperProperties();
    props[0]   = kc;
    props[9]   = uDilMax;
    props.push_back( txCrit );
    const CohesiveFrictionalCrackSubmodels m( props.data(), props.size(), 1 );

    for ( const double uT : { 1e-8, 1e-3, 0.1, 10.0 } ) {
      const double uDil1  = m.uDilMax * std::tanh( m.eta / m.uDilMax * uT );
      double       txLast = -std::numeric_limits< double >::infinity();
      // from deep in the plateau to the contact limit, with penetrations down to round-off
      std::vector< double > uxs;
      for ( int i = 0; i <= 200; i++ )
        uxs.push_back( uDil1 - 4 * std::abs( txCrit ) / m.Kc * std::pow( 10.0, -8.0 * i / 200 ) );
      for ( const double ux : uxs ) {
        const double tx    = m.solveContactTraction( ux, uT );
        const double scale = m.Kc * ( std::abs( ux ) + uDil1 );
        const auto   p     = m.pressureFunction( tx );
        const double th    = std::tanh( m.eta * p.g / m.uDilMax * uT );
        const double r     = tx - m.Kc * ( ux - m.uDilMax * th );
        // the error of t_x: the residual over its slope dr/dt_x = 1 + K_c du_dil/dg dg/dt_x
        const double      error = std::abs( r ) / ( 1.0 + m.Kc * m.eta * uT * ( 1 - th * th ) * p.dG_dTx );
        const std::string where = Marmot::Utilities::MakeString()
                                  << "K_c = " << m.Kc << ", u_T = " << uT << ", u_x = " << ux;
        if ( m.Kc * ux <= txCrit )
          throwExceptionOnFailure( std::abs( tx - m.Kc * ux ) <= 1e-14 * scale,
                                   Marmot::Utilities::MakeString()
                                     << where << ": on the plateau t_x = K_c u_x, but " << tx );
        throwExceptionOnFailure( tx >= std::min( txCrit, m.Kc * ux ) - 1e-14 * scale && tx <= 0,
                                 Marmot::Utilities::MakeString() << where << ": t_x = " << tx << " out of range" );
        throwExceptionOnFailure( error <= 2e-14 * scale,
                                 Marmot::Utilities::MakeString()
                                   << where << ": error " << error << " (residual " << r << ") at t_x = " << tx );
        throwExceptionOnFailure( tx >= txLast - 2e-14 * scale,
                                 Marmot::Utilities::MakeString() << where << ": t_x " << tx << " decreases with u_x" );
        txLast = tx;
      }
    }
  }

  auto props = paperProperties();
  props.push_back( -3.0 );
  const CohesiveFrictionalCrackSubmodels m( props.data(), props.size(), 1 );
  bool                                   thrown = false;
  try {
    m.solveContactTraction( std::numeric_limits< double >::quiet_NaN(), 0.1 );
  }
  catch ( const std::invalid_argument& ) {
    thrown = true;
  }
  throwExceptionOnFailure( thrown, "a non-finite jump must throw" );
}

void testExtrapolation()
{
  const auto props = paperProperties();

  // in the first increment nothing is extrapolated: a jump beyond the strength is elastic
  {
    CohesivePoint p( props );
    const auto&   par = p.model();
    p.moveTo( Eigen::Vector3d( 3 * par.alpha0, 0, 0 ) );
    throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 3 * par.ft, 1e-12 * par.ft ), "elastic first increment" );
    throwExceptionOnFailure( p.getState( "damage" ) > 0, "the implicit damage grows nevertheless" );
  }

  // the damage of an increment is extrapolated linearly in time, alpha_n + dT / dT_n ( alpha_n - alpha_n-1 )
  {
    CohesivePoint p( props );
    const auto&   par   = p.model();
    const double  du    = 0.2 * par.alpha0;
    const double  uxEnd = 0.5 * ( par.alpha0 + par.alphaInt );
    loadTo( p, Eigen::Vector3d( uxEnd, 0, 0 ), int( uxEnd / du ) );
    const double alphaN = p.getState( "alphaMax" ), alphaNm1 = p.getState( "alphaMaxOld" );
    throwExceptionOnFailure( checkIfEqual( alphaN - alphaNm1, uxEnd / int( uxEnd / du ), 1e-9 * du ),
                             "history of the extrapolation" );

    for ( const auto& [dT, dux, what] :
          std::vector< std::tuple< double, double, std::string > >{ { 1.0, 0.0, "same time increment" },
                                                                    { 2.0, 0.0, "double time increment" },
                                                                    { 0.5,
                                                                      -0.1 * par.alpha0,
                                                                      "half time increment, unloading" } } ) {
      auto q          = p.copy( Update::implEx );
      q.dT            = dT;
      const double ux = p.jump( 0 ) + dux;
      q.moveTo( Eigen::Vector3d( ux, 0, 0 ) );
      const double expected = ( 1 - damageOf( par, alphaN + dT * ( alphaN - alphaNm1 ) ) ) * par.Kc * ux;
      throwExceptionOnFailure( checkIfEqual( q.traction( 0 ), expected, 1e-10 * par.ft ),
                               Marmot::Utilities::MakeString()
                                 << what << ": " << q.traction( 0 ) << " vs " << expected );
    }

    // an increment without time extrapolates nothing and keeps the history of the extrapolation
    auto z = p.copy( Update::implEx );
    z.dT   = 0.0;
    z.moveTo( p.jump );
    throwExceptionOnFailure( checkIfEqual( z.traction( 0 ),
                                           ( 1 - damageOf( par, alphaN ) ) * par.Kc * p.jump( 0 ),
                                           1e-10 * par.ft ),
                             "no extrapolation for dT = 0" );
    throwExceptionOnFailure( z.getState( "alphaMaxOld" ) == alphaNm1 && z.getState( "dTOld" ) == 1.0,
                             "history kept for dT = 0" );

    // a tiny remainder increment, as the one that closes a step after the rounding of the summed time increments,
    // does not shift the history: the next increment extrapolates from the last proper one, and does not amplify
    // the round-off of the remainder by dT / dT_remainder
    auto s = p.copy( Update::implEx );
    s.dT   = 1e-9;
    s.moveTo( p.jump + Eigen::Vector3d( 1e-6 * du, 0, 0 ) );
    throwExceptionOnFailure( s.getState( "alphaMaxOld" ) == alphaNm1 && s.getState( "dTOld" ) == 1.0,
                             "history kept for a tiny remainder increment" );
    const double alphaR = s.getState( "alphaMax" );
    throwExceptionOnFailure( alphaR > alphaN, "the implicit history follows the remainder increment" );
    s.dT             = 1.0;
    const double uxR = s.jump( 0 );
    s.moveTo( s.jump );
    const double expected = ( 1 - damageOf( par, alphaR + ( alphaR - alphaNm1 ) ) ) * par.Kc * uxR;
    throwExceptionOnFailure( checkIfEqual( s.traction( 0 ), expected, 1e-10 * par.ft ),
                             Marmot::Utilities::MakeString() << "extrapolation after a tiny remainder increment: "
                                                             << s.traction( 0 ) << " vs " << expected );
  }
}

void testImplExModeIPeakAndEnergy()
{
  // With the jump as time, the extrapolation of a proportional path is exact apart from the kinks of the softening
  // law. The first increment, which is not extrapolated, stays in the elastic range.
  CohesivePoint p( paperProperties() );
  const auto&   par = p.model();

  const int    nSteps = 200000;
  const double uxEnd  = 1.2 * par.alphaUlt;

  double peak = 0, energy = 0;
  p.dT = 0.5 * par.alpha0;
  p.moveTo( Eigen::Vector3d( p.dT, 0, 0 ) );
  energy = 0.5 * p.traction( 0 ) * p.jump( 0 );
  p.dT   = uxEnd / nSteps;
  for ( int i = 1; i <= nSteps; i++ ) {
    const double tOld = p.traction( 0 );
    const double uOld = p.jump( 0 );
    p.moveTo( Eigen::Vector3d( 0.5 * par.alpha0 + uxEnd * i / nSteps, 0, 0 ) );
    peak = std::max( peak, p.traction( 0 ) );
    energy += 0.5 * ( tOld + p.traction( 0 ) ) * ( p.jump( 0 ) - uOld );
  }

  // the area under the bilinear law of Eq. 19, including the (recoverable) elastic branch
  const double area = 0.5 * par.ft * par.alpha0 + par.GF - 0.5 * par.ftInt * ( par.alpha1 - par.alphaInt );

  throwExceptionOnFailure( checkIfEqual( peak, par.ft, 1e-4 * par.ft ),
                           Marmot::Utilities::MakeString() << "peak " << peak );
  throwExceptionOnFailure( checkIfEqual( energy, area, 1e-4 * area ),
                           Marmot::Utilities::MakeString() << "energy " << energy << " vs " << area );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.0, 1e-12 ), "traction after full separation" );
  throwExceptionOnFailure( checkIfEqual( p.getState( "damage" ), 1.0, 1e-14 ), "damage after full separation" );
  throwExceptionOnFailure( checkIfEqual( p.elasticEnergyDensity, 0.0, 1e-14 ), "stored energy after full separation" );
  throwExceptionOnFailure( checkIfEqual( p.dissipation, energy, 1e-12 * energy ),
                           Marmot::Utilities::MakeString()
                             << "dissipation " << p.dissipation << " vs work " << energy );
}

void testImplExSofteningBranchesAndUnloading()
{
  CohesivePoint p( paperProperties() );
  const auto&   par = p.model();

  // on the first branch, the traction follows Eq. 19
  const double ux1 = 0.5 * ( par.alpha0 + par.alphaInt );
  loadTo( p, Eigen::Vector3d( ux1, 0, 0 ), 100 );
  const double q1 = par.ft * ( par.alpha1 - ux1 ) / ( par.alpha1 - par.alpha0 );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), q1, 1e-10 ), "first softening branch" );

  // The first unloading increment still extrapolates the loading: its damage is the one of the next loading
  // increment. The implicit damage is frozen, and from the second unloading increment on the secant is exact.
  const double d          = p.getState( "damage" );
  const double dAlpha     = p.getState( "alphaMax" ) - p.getState( "alphaMaxOld" );
  const double dOvershoot = damageOf( par, ux1 + dAlpha );
  p.moveTo( Eigen::Vector3d( 0.5 * ux1, 0, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), ( 1 - dOvershoot ) * par.Kc * 0.5 * ux1, 1e-10 ),
                           "first unloading increment" );
  throwExceptionOnFailure( checkIfEqual( p.getState( "damage" ), d, 1e-15 ), "damage frozen during unloading" );

  p.moveTo( Eigen::Vector3d( 0.4 * ux1, 0, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), ( 1 - d ) * par.Kc * 0.4 * ux1, 1e-10 ), "secant unloading" );
  throwExceptionOnFailure( checkIfEqual( p.tangent( 0, 0 ), ( 1 - d ) * par.Kc, 1e-6 ), "secant stiffness" );

  // the secant unloading dissipates nothing
  const double dissipation = p.dissipation;
  loadTo( p, Eigen::Vector3d( 0.1 * ux1, 0, 0 ), 5 );
  throwExceptionOnFailure( checkIfEqual( p.dissipation, dissipation, 1e-12 * dissipation ),
                           "unloading must not dissipate" );

  // on the second branch
  const double ux2 = 0.5 * ( par.alphaInt + par.alphaUlt );
  loadTo( p, Eigen::Vector3d( ux2, 0, 0 ), 100 );
  const double q2 = par.ftInt * ( par.alphaUlt - ux2 ) / ( par.alphaUlt - par.alphaInt );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), q2, 1e-10 ), "second softening branch" );
}

void testImplExFrictionalSliding()
{
  auto props = paperProperties();
  props[8]   = 0.0; // no dilatancy

  CohesivePoint p( props );
  const auto&   par = p.model();

  // open the crack completely, then close it under compression and shear it; the faces touch in the last closing
  // increment, and the frozen contact activates in the next one
  loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  const double uxClosed = -1e-7;
  loadTo( p, Eigen::Vector3d( uxClosed, 0, 0 ), 100 );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.0, 1e-9 ), "frozen contact in the closing increment" );
  p.moveTo( p.jump );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), par.Kc * uxClosed, 1e-9 ), "penalty contact" );

  // in steady sliding with constant increments, the extrapolated slip is exact
  loadTo( p, Eigen::Vector3d( uxClosed, 0.3, 0.4 ), 100 );
  const double limit = par.mu * par.Kc * -uxClosed;
  const double tT    = p.traction.tail< 2 >().norm();
  throwExceptionOnFailure( checkIfEqual( tT, limit, 1e-9 ), Marmot::Utilities::MakeString() << "sliding " << tT );
  throwExceptionOnFailure( checkIfEqual( p.traction( 2 ) / p.traction( 1 ), 0.4 / 0.3, 1e-10 ), "sliding direction" );
  throwExceptionOnFailure( checkIfEqual( p.getState( "slip", 1 ) / p.getState( "slip", 0 ), 0.4 / 0.3, 1e-10 ),
                           "slip direction" );

  // the frictional work of steady sliding is dissipated
  const double dissipation = p.dissipation;
  loadTo( p, Eigen::Vector3d( uxClosed, 0.6, 0.8 ), 100 );
  throwExceptionOnFailure( checkIfEqual( p.dissipation - dissipation, limit * 0.5, 1e-9 * limit ),
                           Marmot::Utilities::MakeString()
                             << "frictional dissipation " << p.dissipation - dissipation << " vs " << limit * 0.5 );

  // On reversal, the extrapolated slip of the forward sliding exceeds the trial sliding: the frictional traction
  // vanishes for one increment; then the slip is extrapolated from the implicit (sticking) increment
  p.moveTo( Eigen::Vector3d( uxClosed, 0.6 - 1e-9, 0.8 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction.tail< 2 >().norm(), 0.0, 1e-12 ), "first reversal increment" );
  loadTo( p, Eigen::Vector3d( uxClosed, -0.3, -0.4 ), 100 );
  throwExceptionOnFailure( checkIfEqual( p.traction.tail< 2 >().norm(), limit, 1e-9 ) && p.traction( 1 ) < 0,
                           "reverse sliding" );
}

void testImplExDilatancy()
{
  CohesivePoint p( paperProperties() );
  const auto&   par = p.model();

  loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  loadTo( p, Eigen::Vector3d( 0, 0, 0 ), 100 );

  // the dilatancy of the first shear increment is the one of the converged tangential jump, zero
  p.moveTo( Eigen::Vector3d( 0, 0.02, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.0, 1e-12 ), "dilatancy lags one increment" );

  // with constant increments, the extrapolated tangential jump is exact: the prevented dilatancy produces compression
  const double uy = 2.0;
  loadTo( p, Eigen::Vector3d( 0, uy, 0 ), 99 );
  const double uDil = par.uDilMax * std::tanh( par.eta / par.uDilMax * uy );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), -par.Kc * uDil, 1e-9 * par.Kc * uDil ),
                           "dilatancy pressure" );

  // opened beyond the dilatancy curve (u_T held, so that its extrapolation is exact), the faces separate; with the
  // contact frozen at the last converged increment, they transmit the tension K_c ( u_x - u_dil ) for one increment
  p.moveTo( Eigen::Vector3d( 0, uy, 0 ) );
  const double uxOpen = uDil + 1e-9;
  p.moveTo( Eigen::Vector3d( uxOpen, uy, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), par.Kc * ( uxOpen - uDil ), 1e-9 ),
                           Marmot::Utilities::MakeString() << "frozen contact on opening: " << p.traction( 0 ) );
  throwExceptionOnFailure( p.getState( "contact" ) == 0.0, "contact released after opening" );
  p.moveTo( Eigen::Vector3d( uxOpen, uy, 0 ) );
  throwExceptionOnFailure( checkIfEqual( p.traction.norm(), 0.0, 1e-9 ), "free dilatant opening" );
}

void testImplExDilatancyPressureCutoff()
{
  auto props = paperProperties();
  props.push_back( -1.0 ); // t_x,crit

  CohesivePoint p( props );
  const auto&   par = p.model();

  loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 2000 );
  loadTo( p, Eigen::Vector3d( -1e-6, 0, 0 ), 100 ); // normal traction -34 < t_x,crit
  loadTo( p, Eigen::Vector3d( -1e-6, 0.01, 0 ), 10 );

  // t_x = -34 < t_x,crit: g = 0, no dilatancy
  throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), -par.Kc * 1e-6, 1e-9 ), "dilatancy pressure cutoff" );

  // between t_x,crit and 0, g( t_x ) is evaluated with the normal traction it produces: for u_x = 0 and a small
  // (extrapolated) tangential jump, t_x = -a g = -a ( 1 + t_x ) with a = K_c u_dil(g = 1)
  CohesivePoint fresh( props );
  const double  uT = 0.5 / ( par.Kc * par.eta );
  fresh.moveTo( Eigen::Vector3d( 0, 0.5 * uT, 0 ) );
  fresh.moveTo( Eigen::Vector3d( 0, uT, 0 ) );
  const double a = par.Kc * par.uDilMax * std::tanh( par.eta / par.uDilMax * uT );
  throwExceptionOnFailure( checkIfEqual( fresh.traction( 0 ), -a / ( 1 + a ), 1e-9 ),
                           Marmot::Utilities::MakeString()
                             << "implicit pressure function: " << fresh.traction( 0 ) << " vs " << -a / ( 1 + a ) );
  checkTangent( fresh, Eigen::Vector3d( -1e-9, 1.1 * uT, 0.2 * uT ), "implicit pressure function" );
}

void testImplExTangents()
{
  {
    CohesivePoint p( paperProperties() );
    const auto&   par = p.model();
    checkTangent( p, Eigen::Vector3d( 0.3, 0.2, -0.1 ) * par.alpha0, "elastic" );

    loadTo( p, Eigen::Vector3d( 2.0, 1.0, 0.5 ) * par.alpha0, 50 );
    checkTangent( p, Eigen::Vector3d( 2.0, 1.0, 0.5 ) * par.alpha0 * 1.01, "softening, mixed mode" );
    checkTangent( p, Eigen::Vector3d( 2.0, 1.0, 0.5 ) * par.alpha0 * 0.9, "unloading" );

    loadTo( p, Eigen::Vector3d( 0.3 * par.alphaUlt, 0.01, 0 ), 50 );
    checkTangent( p, Eigen::Vector3d( 0.31 * par.alphaUlt, 0.012, 0.001 ), "second branch, mixed mode" );
    p.dT = 3.0;
    checkTangent( p, Eigen::Vector3d( 0.31 * par.alphaUlt, 0.012, 0.001 ), "second branch, larger time increment" );
  }
  {
    CohesivePoint p( paperProperties() );
    const auto&   par = p.model();
    loadTo( p, Eigen::Vector3d( 0.2 * par.alphaUlt, 0, 0 ), 50 );
    loadTo( p, Eigen::Vector3d( -1e-6, 0, 0 ), 20 );
    checkTangent( p, Eigen::Vector3d( -1e-6, 1e-7, 0 ), "partially damaged, contact stick" );
    loadTo( p, Eigen::Vector3d( -1e-6, 0.05, 0.02 ), 20 );
    checkTangent( p, Eigen::Vector3d( -1e-6, 0.053, 0.021 ), "partially damaged, contact slip with dilatancy" );
  }
  {
    CohesivePoint p( paperProperties() );
    const auto&   par = p.model();
    loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 100 );
    checkTangent( p, Eigen::Vector3d( 1.2 * par.alphaUlt, 0.01, 0 ), "completely separated, open" );
    loadTo( p, Eigen::Vector3d( -1e-6, 0, 0 ), 20 );
    loadTo( p, Eigen::Vector3d( -1e-6, 0.3, -0.1 ), 20 );
    checkTangent( p, Eigen::Vector3d( -1e-6, 0.315, -0.105 ), "completely separated, contact slip" );
    checkTangent( p, Eigen::Vector3d( -1e-6, 0.3, -0.1 ), "completely separated, extrapolated slip beyond trial" );
  }
  {
    // pressure cutoff with t_x,crit < t_x < 0, where g( t_x ) is implicit
    auto props = paperProperties();
    props[9]   = 0.1; // saturating dilatancy
    props.push_back( -3.0 );
    CohesivePoint p( props );
    const auto&   par = p.model();
    loadTo( p, Eigen::Vector3d( 0.2 * par.alphaUlt, 0, 0 ), 50 );
    loadTo( p, Eigen::Vector3d( 0.0, 0, 0 ), 20 );
    loadTo( p, Eigen::Vector3d( 0.0, 2e-6, 1e-6 ), 2 );
    checkTangent( p, Eigen::Vector3d( 0.0, 3e-6, 1.5e-6 ), "partially damaged, implicit pressure function, stick" );
    loadTo( p, Eigen::Vector3d( 0.0, 0.05, 0.02 ), 20 );
    checkTangent( p, Eigen::Vector3d( 0.0, 0.0525, 0.021 ), "partially damaged, implicit pressure function, slip" );
  }
}

/// Impl-Ex freezes the contact indicator at the last converged increment; the explicit update follows the current
/// jump.
void testContact()
{
  const double uxClosed = -1e-6, uxOpen = 1e-6;

  for ( const auto update : { Update::implEx, Update::explicitDynamics } ) {
    const bool    frozen = update == Update::implEx;
    CohesivePoint p( paperProperties(), update );
    const auto&   par  = p.model();
    const auto    name = frozen ? "Impl-Ex" : "explicit";

    // a completely separated crack, open; then closed in one increment and opened again in one increment
    loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 100 );
    throwExceptionOnFailure( p.getState( "contact" ) == 0.0, "open after separation" );

    // closing: frozen, the contact activates one increment later, the faces penetrate without resistance
    p.moveTo( Eigen::Vector3d( uxClosed, 0, 0 ) );
    throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), frozen ? 0.0 : par.Kc * uxClosed, 1e-9 ),
                             Marmot::Utilities::MakeString() << name << ", closing increment: " << p.traction( 0 ) );
    throwExceptionOnFailure( p.getState( "contact" ) == 1.0, "contact stored after closing" );
    p.moveTo( Eigen::Vector3d( uxClosed, 0, 0 ) );
    throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), par.Kc * uxClosed, 1e-9 ),
                             Marmot::Utilities::MakeString() << name << ", in contact: " << p.traction( 0 ) );

    // opening: frozen, the contact transmits tension for one increment
    p.moveTo( Eigen::Vector3d( uxOpen, 0, 0 ) );
    throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), frozen ? par.Kc * uxOpen : 0.0, 1e-9 ),
                             Marmot::Utilities::MakeString() << name << ", opening increment: " << p.traction( 0 ) );
    throwExceptionOnFailure( p.getState( "contact" ) == 0.0, "contact released after opening" );
    p.moveTo( Eigen::Vector3d( uxOpen, 0, 0 ) );
    throwExceptionOnFailure( checkIfEqual( p.traction( 0 ), 0.0, 1e-9 ),
                             Marmot::Utilities::MakeString() << name << ", open: " << p.traction( 0 ) );
  }
}

/**
 * The maximum stiffness bounds the eigenvalues of the Impl-Ex tangent in all states, including contact with sliding
 * and strong dilatancy.
 */
void testMaximumStiffness()
{
  auto strongDilatancy = paperProperties();
  strongDilatancy[7]   = 0.9; // mu
  strongDilatancy[8]   = 0.9; // eta

  for ( const auto& props : { paperProperties(), strongDilatancy } ) {
    CohesivePoint p( props );
    const auto&   par = p.model();
    // the bound does not depend on the state: from the virgin state, it must also hold for all later ones
    const Marmot::Materials::MaterialCohesive::Response virgin{ p.traction,
                                                                p.stateVars.data(),
                                                                p.elasticEnergyDensity,
                                                                p.dissipation };
    const double                                        maxStiffness = par.getMaximumStiffness( virgin, p.jump );

    throwExceptionOnFailure( checkIfEqual( maxStiffness,
                                           par.Kc * ( std::max( 1.0, par.beta * par.beta / par.kappa ) +
                                                      std::sqrt( par.mu * par.eta ) ),
                                           1e-12 * par.Kc ),
                             "maximum stiffness" );

    // a trial jump from the current history: the eigenvalues of the tangent
    const auto check = [&]( const Eigen::Vector3d& trialJump, const std::string& where ) {
      Marmot::Materials::MaterialCohesive::Tangents tangents;
      p.update( trialJump, tangents, false );
      const double spectralRadius = tangents.dTraction_dJump.eigenvalues().cwiseAbs().maxCoeff();
      throwExceptionOnFailure( spectralRadius <= maxStiffness,
                               Marmot::Utilities::MakeString()
                                 << where << ": an eigenvalue of the tangent, " << spectralRadius
                                 << ", exceeds the maximum stiffness " << maxStiffness );
    };

    check( Eigen::Vector3d( 0.3, 0.2, -0.1 ) * par.alpha0, "elastic" );
    loadTo( p, Eigen::Vector3d( 2.0, 1.0, 0.5 ) * par.alpha0, 50 );
    check( Eigen::Vector3d( 2.0, 1.0, 0.5 ) * par.alpha0 * 1.01, "softening" );
    loadTo( p, Eigen::Vector3d( 0.2 * par.alphaUlt, 0, 0 ), 50 );
    loadTo( p, Eigen::Vector3d( -1e-6, 0, 0 ), 20 );
    check( Eigen::Vector3d( -1e-6, 1e-7, 0 ), "partially damaged, contact stick" );
    check( Eigen::Vector3d( -1e-6, 0.05, 0.02 ), "partially damaged, contact slip with dilatancy" );
    loadTo( p, Eigen::Vector3d( 1.1 * par.alphaUlt, 0, 0 ), 100 );
    loadTo( p, Eigen::Vector3d( -1e-6, 0, 0 ), 20 );
    check( Eigen::Vector3d( -1e-6, 0.3, -0.1 ), "completely separated, contact slip" );
    check( Eigen::Vector3d( -1e-6, 1e-8, 0 ), "completely separated, contact stick" );
    loadTo( p, Eigen::Vector3d( -1e-3, 1e-2, 0 ), 20 );
    check( Eigen::Vector3d( -1e-3, 1.2e-2, 1e-3 ), "completely separated, deep contact, slip with dilatancy" );
  }
}

void testInvalidParameters()
{
  auto expectThrow = []( std::vector< double > props, const std::string& what ) {
    bool thrown = false;
    try {
      CohesiveFrictionalCrack( props.data(), props.size(), 1 );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown, "no exception for " + what );
  };

  auto props = paperProperties();
  expectThrow( { props.begin(), props.begin() + 9 }, "too few properties" );
  props[4] = 3.0;
  expectThrow( props, "f_t,int > f_t" );
  props    = paperProperties();
  props[3] = 0.2;
  expectThrow( props, "G_f > G_F" );
  props = paperProperties();
  props.insert( props.end(), { 0.0, 0.0, 1.7e-9 } );
  expectThrow( props, "the former layout with the tangent method at [11] and the mass per area at [12]" );
  props = paperProperties();
  props.insert( props.end(), { 0.0, -1e-9 } );
  expectThrow( props, "negative mass per area" );
  props = paperProperties();
  props.insert( props.end(), { 0.0, std::nan( "" ) } );
  expectThrow( props, "NaN mass per area" );
}

/// The mass per area is the optional last parameter; as the density of a bulk material, asking for it throws if it
/// was not given.
void testMassPerArea()
{
  auto props  = paperProperties();
  bool thrown = false;
  try {
    CohesiveFrictionalCrack( props.data(), props.size(), 1 ).getMassPerArea( nullptr );
  }
  catch ( const std::runtime_error& ) {
    thrown = true;
  }
  throwExceptionOnFailure( thrown, "asking for a mass per area that was not given must throw" );

  const double massPerArea = 1.7e-9;
  props.insert( props.end(), { 0.0, massPerArea } );
  throwExceptionOnFailure( CohesiveFrictionalCrack( props.data(), props.size(), 1 ).getMassPerArea( nullptr ) ==
                             massPerArea,
                           "the mass per area is materialProperties[11]" );
}

int main()
{
  executeTestsAndCollectExceptions( { testExplicitModeIPeakAndEnergy,
                                      testExplicitEnergies,
                                      testExplicitSofteningBranchesAndUnloading,
                                      testExplicitFrictionalSliding,
                                      testExplicitDilatancy,
                                      testExplicitDilatancyPressureCutoff,
                                      testExplicitDoesNotExtrapolate,
                                      testImplExAndExplicitStates,
                                      testShearStrength,
                                      testEquivalentOpening,
                                      testSubmodelDerivatives,
                                      testContactTraction,
                                      testExtrapolation,
                                      testImplExModeIPeakAndEnergy,
                                      testImplExSofteningBranchesAndUnloading,
                                      testImplExFrictionalSliding,
                                      testImplExDilatancy,
                                      testImplExDilatancyPressureCutoff,
                                      testImplExTangents,
                                      testContact,
                                      testMaximumStiffness,
                                      testInvalidParameters,
                                      testMassPerArea } );
  return 0;
}
