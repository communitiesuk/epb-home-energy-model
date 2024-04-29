use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::space_heat_demand::building_element::projected_height;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::external_conditions::ExternalConditions;
use crate::input::OnSiteGenerationVentilationStrategy;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

/// This module contains objects that represent photovoltaic systems.

/// system performance factor lookup
/// informative values from table C.4 Annex C BS EN 15316-4-3:2017
/// note from 6.2.4.7.2 rear surface free - is if PV system is not integrated.
/// Assume this means NOT integrated (BIPV)  or attached (BAPV)
/// Increased by 0.85/0.8 based on median quoted performance ratio from
/// "Performance of Distributed PV in the UK: A Statistical Analysis of
/// Over 7000 systems" conference paper from 31st European Photovoltaic
/// Solar Energy Conference and Exhibition, September 2015, assuming this
/// applies to moderately ventilated case.
// Note: BS EN 15316-4-3:2017 section 6.2.4.7.2 states that the performance
// factor for "rear surface free" should be 1.0. However, this would
// seem to imply that there are no inverter or other system losses,
// despite the fact that section 6.2.4.7.5 states that this factor
// accounts for these losses. Also, no factor for "rear surface free"
// has been given in Table C.4. Therefore, it was decided to use the
// same factor for "rear surface free" as for "strongly or forced
// ventilated".
const F_PERF_LOOKUP_UNVENTILATED: f64 = 0.81;
const F_PERF_LOOKUP_MODERATELY_VENTILATED: f64 = 0.85;
const F_PERF_LOOKUP_STRONGLY_OR_FORCED_VENTILATED: f64 = 0.87;
const F_PERF_LOOKUP_REAR_SURFACE_FREE: f64 = 0.87;

pub struct PhotovoltaicSystem {
    peak_power: f64,
    f_perf: f64,
    pitch: f64,
    orientation: f64,
    base_height: f64,
    width: f64,
    projected_height: f64,
    external_conditions: Arc<ExternalConditions>,
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
}

impl PhotovoltaicSystem {
    /// Construct a PhotovoltaicSystem object
    ///
    /// Arguments:
    /// * `peak_power` - Peak power in kW; represents the electrical power of a photovoltaic
    ///                  system with a given area and a for a solar irradiance of 1 kW/m2
    ///                  on this surface (at 25 degrees)
    ///                  TODO (from Python) - Could add other options at a later stage.
    ///                  Standard has alternative method when peak power is not available
    ///                  (input type of PV module and Area instead when peak power unknown)
    /// * `ventilation_strategy` - ventilation strategy of the PV system.
    ///                            This will be used to determine the system performance factor
    ///                            based on a lookup table
    /// * `pitch` - is the tilt angle (inclination) of the PV panel from horizontal,
    ///             measured upwards facing, 0 to 90, in degrees.
    ///             0=horizontal surface, 90=vertical surface.
    ///             Needed to calculate solar irradiation at the panel surface.
    /// * `orientation` - is the orientation angle of the inclined surface, expressed as the
    ///                   geographical azimuth angle of the horizontal projection of the inclined
    ///                   surface normal, -180 to 180, in degrees;
    ///                   Assumed N 180 or -180, E 90, S 0, W -90
    ///                   TODO (from Python) - PV standard refers to angle as between 0 to 360?
    ///                   Needed to calculate solar irradiation at the panel surface.
    /// * `base_height` - is the distance between the ground and the lowest edge of the PV panel, in m
    /// * `height` - is the height of the PV panel, in m
    /// * `width` - is the width of the PV panel, in m
    /// * `external_conditions` - reference to ExternalConditions object
    /// * `energy_supply_connection` - an EnergySupplyConnection value
    /// * `simulation_timestep` - reference to step length of a SimulationTime object in the context
    pub fn new(
        peak_power: f64,
        ventilation_strategy: OnSiteGenerationVentilationStrategy,
        pitch: f64,
        orientation: f64,
        base_height: f64,
        height: f64,
        width: f64,
        external_conditions: Arc<ExternalConditions>,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
    ) -> Self {
        Self {
            peak_power,
            f_perf: match ventilation_strategy {
                OnSiteGenerationVentilationStrategy::Unventilated => F_PERF_LOOKUP_UNVENTILATED,
                OnSiteGenerationVentilationStrategy::ModeratelyVentilated => {
                    F_PERF_LOOKUP_MODERATELY_VENTILATED
                }
                OnSiteGenerationVentilationStrategy::StronglyOrForcedVentilated => {
                    F_PERF_LOOKUP_STRONGLY_OR_FORCED_VENTILATED
                }
                OnSiteGenerationVentilationStrategy::RearSurfaceFree => {
                    F_PERF_LOOKUP_REAR_SURFACE_FREE
                }
            },
            pitch,
            orientation,
            base_height,
            width,
            projected_height: projected_height(pitch, height),
            external_conditions,
            energy_supply_connection,
            simulation_timestep,
        }
    }

    /// Produce electrical energy (in kWh) from the PV system
    /// according to BS EN 15316-4-3:2017
    pub fn produce_energy(&self, simulation_time_iteration: SimulationTimeIteration) {
        // solar irradiance in W/m2
        let (i_sol_dir, i_sol_dif, _, _) = self
            .external_conditions
            .calculated_direct_diffuse_total_irradiance(
                self.pitch,
                self.orientation,
                false,
                &simulation_time_iteration,
            );
        // shading factors
        let (f_sh_dir, f_sh_dif) = self.shading_factors_direct_diffuse(simulation_time_iteration);
        // solar irradiation in kWh/m2
        let solar_irradiation = (i_sol_dir * f_sh_dir + i_sol_dif * f_sh_dif)
            * self.simulation_timestep
            / WATTS_PER_KILOWATT as f64;
        // reference solar irradiance kW/m2
        let ref_solar_irradiance = 1.;

        // CALCULATION
        // E.el.pv.out.h = E.sol.pv.h * P.pk * f.perf / I.ref
        // energy_produced = solar_irradiation * peak_power * system_performance_factor
        //                     / reference_solar_irradiance
        // energy produced in kWh
        let energy_produced =
            solar_irradiation * self.peak_power * self.f_perf / ref_solar_irradiance;

        self.energy_supply_connection
            .supply_energy(energy_produced, simulation_time_iteration.index)
            .unwrap();
    }

    fn shading_factors_direct_diffuse(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        self.external_conditions
            .shading_reduction_factor_direct_diffuse(
                self.base_height,
                self.projected_height,
                self.width,
                self.pitch,
                self.orientation,
                vec![],
                simulation_time_iteration,
            )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use parking_lot::Mutex;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    pub fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
            51.42,
            -0.75,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            vec![
                ShadingSegment {
                    number: 1,
                    start: 180,
                    end: 135,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 2,
                    start: 135,
                    end: 90,
                    objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Overhang,
                        height: 2.2,
                        distance: 6.,
                    }]),
                },
                ShadingSegment {
                    number: 3,
                    start: 90,
                    end: 45,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 4,
                    start: 45,
                    end: 0,
                    objects: Some(vec![
                        ShadingObject {
                            object_type: ShadingObjectType::Obstacle,
                            height: 40.,
                            distance: 4.,
                        },
                        ShadingObject {
                            object_type: ShadingObjectType::Overhang,
                            height: 3.,
                            distance: 7.,
                        },
                    ]),
                },
                ShadingSegment {
                    number: 5,
                    start: 0,
                    end: -45,
                    objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Obstacle,
                        height: 3.,
                        distance: 8.,
                    }]),
                },
                ShadingSegment {
                    number: 6,
                    start: -45,
                    end: -90,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 7,
                    start: -90,
                    end: -135,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 8,
                    start: -135,
                    end: -180,
                    ..Default::default()
                },
            ],
        )
    }

    #[fixture]
    pub fn pv(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) -> (PhotovoltaicSystem, Arc<Mutex<EnergySupply>>) {
        let energy_supply = Arc::new(Mutex::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
        )));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "pv generation").unwrap();
        let pv = PhotovoltaicSystem::new(
            2.5,
            OnSiteGenerationVentilationStrategy::ModeratelyVentilated,
            30.,
            0.,
            10.,
            2.,
            3.,
            Arc::new(external_conditions),
            energy_supply_conn,
            simulation_time.step,
        );
        (pv, energy_supply)
    }

    #[rstest]
    pub fn test_produce_energy(
        pv: (PhotovoltaicSystem, Arc<Mutex<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, energy_supply) = pv;
        let expected_generation_results = [
            -0.019039734375,
            -0.04586317375,
            -0.072234285625,
            -0.08283324375,
            -0.089484801875,
            -0.074557144375,
            -0.05198861375,
            -0.040747119375,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            pv.produce_energy(t_it);
            assert_eq!(
                round_by_precision(
                    energy_supply.lock().results_by_end_user()["pv generation"][t_idx],
                    1e6
                ),
                round_by_precision(expected_generation_results[t_idx], 1e6),
                "incorrect electricity produced from pv returned"
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }
}
