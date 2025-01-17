/// This module contains objects that represent photovoltaic systems.
use crate::compare_floats::min_of_2;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::space_heat_demand::building_element::projected_height;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::external_conditions::{ExternalConditions, WindowShadingObject};
use crate::input::OnSiteGenerationVentilationStrategy;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

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

#[derive(Debug)]
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
    shading: Vec<WindowShadingObject>,
    inverter_peak_power: f64,
    inverter_is_inside: bool,
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
    /// * `shading`
    /// * `inverter_peak_power` - Peak power in kW; represents the peak electrical power input to the inverter
    /// * `inverter_is_inside` - tells us that the inverter is considered inside the building
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
        shading: Vec<WindowShadingObject>,
        inverter_peak_power: f64,
        inverter_is_inside: bool,
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
            shading,
            inverter_peak_power,
            inverter_is_inside,
        }
    }

    /// Return whether this unit is considered inside the building or not
    pub fn inverter_is_inside(&self) -> bool {
        self.inverter_is_inside
    }

    /// Produce electrical energy (in kWh) from the PV system
    /// according to BS EN 15316-4-3:2017
    pub fn produce_energy(&self, simulation_time_iteration: SimulationTimeIteration) -> (f64, f64) {
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

        // energy input in kWh; now need to calculate total energy produce taking into account inverter efficiency
        let energy_input =
            solar_irradiation * self.peak_power * self.f_perf / 0.92 / ref_solar_irradiance;
        // f_perf is divided by 0.92 to avoid double-applying the inverter efficiency,
        // which is applied separately below via 'inverter_dc_ac_efficiency', since
        // inverter efficiency was inherently included in the factors taken
        // from BS EN 15316-4-3:2017.

        // power output from PV panel in kW used to calculate ratio for efficiency loss of inverters from DC to AC
        let power_input_inverter = energy_input / self.simulation_timestep;

        // Calculate Ratio of Rated Power
        let ratio_of_rated_output =
            min_of_2(power_input_inverter, self.inverter_peak_power) / self.inverter_peak_power;

        // Using Ratio of Rated Power, calculate Inverter DC to AC efficiency
        // equation was estimated based on graph from
        // https://www.researchgate.net/publication/260286647_Performance_of_PV_inverters figure 9
        let inverter_dc_ac_efficiency = if ratio_of_rated_output == 0. {
            0.
        } else {
            0.92 * (4.67375 * ratio_of_rated_output).tanh().powf(0.137951)
        };

        // Calculate energy produced output taking into account peak power of inverter + array
        // and inverter DC to AC efficiency
        let energy_produced = min_of_2(
            energy_input,
            self.inverter_peak_power * self.simulation_timestep,
        ) * inverter_dc_ac_efficiency;

        // Add energy produced to the applicable energy supply connection (this will reduce demand)
        self.energy_supply_connection
            .supply_energy(energy_produced, simulation_time_iteration.index)
            .unwrap();

        let energy_lost = energy_input - energy_produced;

        (energy_produced, energy_lost)
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
                &self.shading,
                simulation_time_iteration,
            )
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
        WindowShadingObjectType,
    };
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
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
            vec![220., 230., 240., 250., 260., 270., 270., 280.],
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
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            vec![
                ShadingSegment {
                    number: 1,
                    start: 180.,
                    end: 135.,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 2,
                    start: 135.,
                    end: 90.,
                    shading_objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Overhang,
                        height: 2.2,
                        distance: 6.,
                    }]),
                },
                ShadingSegment {
                    number: 3,
                    start: 90.,
                    end: 45.,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 4,
                    start: 45.,
                    end: 0.,
                    shading_objects: Some(vec![
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
                    start: 0.,
                    end: -45.,
                    shading_objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Obstacle,
                        height: 3.,
                        distance: 8.,
                    }]),
                },
                ShadingSegment {
                    number: 6,
                    start: -45.,
                    end: -90.,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 7,
                    start: -90.,
                    end: -135.,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 8,
                    start: -135.,
                    end: -180.,
                    ..Default::default()
                },
            ],
        )
    }

    #[fixture]
    pub fn pv(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) -> (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "pv generation without shading")
                .unwrap();
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
            vec![],
            2.5,
            false,
        );
        (pv, energy_supply)
    }

    #[fixture]
    pub fn pv_with_shading(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) -> (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "pv generation with shading").unwrap();
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
            vec![
                WindowShadingObject {
                    object_type: WindowShadingObjectType::Overhang,
                    depth: 0.5,
                    distance: 0.5,
                },
                WindowShadingObject {
                    object_type: WindowShadingObjectType::SideFinLeft,
                    depth: 0.25,
                    distance: 0.1,
                },
                WindowShadingObject {
                    object_type: WindowShadingObjectType::SideFinRight,
                    depth: 0.25,
                    distance: 0.1,
                },
            ],
            2.5,
            true,
        );
        (pv, energy_supply)
    }

    #[rstest]
    pub fn test_is_inside(
        pv: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        pv_with_shading: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
    ) {
        let (pv, _) = pv;
        let (pv_with_shading, _) = pv_with_shading;
        assert!(!pv.inverter_is_inside());
        assert!(pv_with_shading.inverter_is_inside());
    }

    #[rstest]
    pub fn test_produce_energy(
        pv: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, energy_supply) = pv;
        let expected_generation_results = [
            -0.012155950159829848,
            -0.033046009462695744,
            -0.05538065916905987,
            -0.06469759901250412,
            -0.07062640777214597,
            -0.057408001837834045,
            -0.038108732702294035,
            -0.028886972604485958,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            pv.produce_energy(t_it);
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation without shading"][t_idx],
                expected_generation_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    pub fn test_produce_energy_with_shading(
        pv_with_shading: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, energy_supply) = pv_with_shading;
        let expected_generation_results = [
            -0.006675561797598833,
            -0.01815144206200631,
            -0.030430978364590064,
            -0.035557634328444936,
            -0.03882148986837187,
            -0.03154628950205703,
            -0.02093382593171027,
            -0.01916405808037656,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            pv.produce_energy(t_it);
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation with shading"][t_idx],
                expected_generation_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    pub fn test_energy_produced_and_energy_lost(
        pv: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, _) = pv;
        let expected_energy_produced = [
            0.012155950159829848,
            0.033046009462695744,
            0.05538065916905987,
            0.06469759901250412,
            0.07062640777214597,
            0.057408001837834045,
            0.038108732702294035,
            0.028886972604485958,
        ];
        let expected_energy_lost = [
            0.008539417678140976,
            0.01680527481710961,
            0.02313487715438004,
            0.02533854258544363,
            0.02663969095724275,
            0.023632380772525823,
            0.018400631461166272,
            0.015403377364203018,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let (energy_produced, energy_lost) = pv.produce_energy(t_it);
            assert_relative_eq!(
                energy_produced,
                expected_energy_produced[t_idx],
                max_relative = 1e-6
            );
            assert_relative_eq!(
                energy_lost,
                expected_energy_lost[t_idx],
                max_relative = 1e-6
            );
        }
    }
}
