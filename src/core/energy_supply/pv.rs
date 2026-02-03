use crate::core::energy_supply::inverter::Inverter;
use crate::core::energy_supply::on_site_generation_base::OnSiteGeneration;
use crate::core::space_heat_demand::building_element::projected_height;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::external_conditions::{
    CalculatedDirectDiffuseTotalIrradiance, ExternalConditions, WindowShadingObject,
};
use crate::input::{InverterType, PhotovoltaicVentilationStrategy};
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

const F_PERF_LOOKUP_UNVENTILATED: f64 = 0.81;
const F_PERF_LOOKUP_MODERATELY_VENTILATED: f64 = 0.85;
const F_PERF_LOOKUP_STRONGLY_OR_FORCED_VENTILATED: f64 = 0.87;
const F_PERF_LOOKUP_REAR_SURFACE_FREE: f64 = 0.87;

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
#[derive(Debug, Clone)]
pub(crate) struct PhotovoltaicPanel {
    peak_power: f64,
    f_perf: f64,
    pitch: f64,
    orientation: f64,
    base_height: f64,
    width: f64,
    projected_height: f64,
    simulation_timestep: f64,
    shading: Vec<WindowShadingObject>,
}

impl PhotovoltaicPanel {
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
    /// * `shading` - TODO (from Python) could add at a later date.
    pub(crate) fn new(
        peak_power: f64,
        ventilation_strategy: PhotovoltaicVentilationStrategy,
        pitch: f64,
        orientation: f64,
        base_height: f64,
        height: f64,
        width: f64,
        simulation_timestep: f64,
        shading: Vec<WindowShadingObject>,
    ) -> Self {
        Self {
            peak_power,
            f_perf: match ventilation_strategy {
                PhotovoltaicVentilationStrategy::Unventilated => F_PERF_LOOKUP_UNVENTILATED,
                PhotovoltaicVentilationStrategy::ModeratelyVentilated => {
                    F_PERF_LOOKUP_MODERATELY_VENTILATED
                }
                PhotovoltaicVentilationStrategy::StronglyOrForcedVentilated => {
                    F_PERF_LOOKUP_STRONGLY_OR_FORCED_VENTILATED
                }
                PhotovoltaicVentilationStrategy::RearSurfaceFree => F_PERF_LOOKUP_REAR_SURFACE_FREE,
            },
            pitch,
            orientation,
            base_height,
            projected_height: projected_height(pitch, height),
            width,
            simulation_timestep,
            shading,
        }
    }

    /// return calculated shading factor
    fn shading_factors_direct(
        &self,
        external_conditions: &Arc<ExternalConditions>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(external_conditions
            .shading_reduction_factor_direct_diffuse(
                self.base_height,
                self.projected_height,
                self.width,
                self.pitch,
                self.orientation,
                &self.shading,
                simulation_time_iteration,
            )?
            .0)
    }

    /// return calculated shading factor
    fn shading_factors_diffuse(
        &self,
        external_conditions: &Arc<ExternalConditions>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(external_conditions
            .shading_reduction_factor_direct_diffuse(
                self.base_height,
                self.projected_height,
                self.width,
                self.pitch,
                self.orientation,
                &self.shading,
                simulation_time_iteration,
            )?
            .1)
    }

    /// Produce electrical energy (in kWh) from the PV system
    /// according to BS EN 15316-4-3:2017
    pub fn produce_energy(
        &self,
        photovoltaic_system: PhotovoltaicSystem,
        external_conditions: &Arc<ExternalConditions>,
        simulation_time_iteration: SimulationTimeIteration,
        f_sh_dir: f64,
    ) -> anyhow::Result<f64> {
        // solar irradiance in W/m2
        let CalculatedDirectDiffuseTotalIrradiance(i_sol_dir, i_sol_dif, _, _) =
            external_conditions.calculated_direct_diffuse_total_irradiance(
                self.pitch,
                self.orientation,
                false,
                &simulation_time_iteration,
            );

        // diffuse shading factor
        let f_sh_dif =
            self.shading_factors_diffuse(external_conditions, simulation_time_iteration)?;

        // Calculate the impact of direct shading on the panel/inverters ability to output energy,
        // i.e. where the shadow of an obstacle falls on the PV panel.
        // There is a lower impact if module level electronics ('optimised inverter') selected.
        // This factor is then applied in the energy_produced calculation later.
        let inv_shad_inefficiency = if f_sh_dir < 1. {
            photovoltaic_system.inverter_efficiency_lookup(f_sh_dir)
        } else {
            1.
        };

        // solar irradiation in kWh/m2
        let solar_irradiation = (i_sol_dir * f_sh_dir + i_sol_dif * f_sh_dif)
            * simulation_time_iteration.timestep
            / WATTS_PER_KILOWATT as f64;

        // reference solar irradiance kW/m2
        let ref_solar_irradiance = 1.;

        // CALCULATION
        // E.el.pv.out.h = E.sol.pv.h * P.pk * f.perf / I. ref
        // energy_produced = solar_irradiation * peak_power * system_performance_factor
        // / reference_solar_irradiance
        //  energy input in kWh; now need to calculate total energy produce taking into account inverter efficiency
        let energy_input =
            solar_irradiation * self.peak_power * self.f_perf * inv_shad_inefficiency
                / 0.972
                / ref_solar_irradiance;
        // (from Python, though out of date comment) f_perf is divided by 0.92 to avoid double-applying the inverter efficiency,
        // which is applied separately below via 'inverter_dc_ac_efficiency', since
        // inverter efficiency was inherently included in the factors taken
        // from BS EN 15316-4-3:2017.

        Ok(energy_input)
    }
}

/// An object to represent a photovoltaic system
#[derive(Clone, Debug)]
pub(crate) struct PhotovoltaicSystem {
    external_conditions: Arc<ExternalConditions>,
    // simulation_timestep: f64,
    panels: Vec<PhotovoltaicPanel>,
    inverter: Inverter,
}
impl PhotovoltaicSystem {
    /// Construct a PhotovoltaicSystem object
    /// Arguments:
    /// * `ext_cond`               -- reference to ExternalConditions object
    /// * `simulation_time`        -- reference to SimulationTime object
    /// * `panels`                 -- a list of PhotovoltaicSystem objects
    /// * `inverter`               -- reference to Inverter object
    pub(crate) fn new(
        external_conditions: Arc<ExternalConditions>,
        panels: Vec<PhotovoltaicPanel>,
        inverter: Inverter,
    ) -> Self {
        Self {
            external_conditions,
            panels,
            inverter,
        }
    }

    ///    Inverter efficiency reduction as a function of shaded proportion and type of inverter is based on data from
    ///    'Partial Shade Evaluation of Distributed Power Electronics for Photovoltaic Systems', Chris Deline et al,
    ///    figure 7, accessed at: https://www.nrel.gov/docs/fy12osti/54039.pdf.
    ///
    ///    The equations given in figure 7 were used to calculate the normalised power output for a range of different
    ///    levels of panel covering, from 0 to 100% for the two inverter types represented (string and micro [aka optimised]).
    ///
    ///    The measured 37% transmission rate of the covering material was applied to determine the reduction in incident
    ///    radiation at each level of coverage:
    ///        % reduction in radiation reaching panel = percent of panel covered * (1 - 0.37)
    ///    (where 1 = all radiation blocked, 0 = no reduction)
    ///
    ///    The level of shading was converted to a shading factor of the form used in HEM
    ///    (where 1 = no shading, 0 = complete shading):
    ///        f_sh_dir = (1 - % reduction in radiation reaching panel)
    ///
    ///    It was assumed that efficiency reduction being calculated is only related to the shading of *direct* radiation
    ///    (not diffuse) on the basis that the impact is caused by part of the panel being shaded, while the rest is not.
    ///    Diffuse shading would affect the whole panel approximately equally, so should not affect inverter efficiency
    ///    in the same way. (Note that the reduction in output associated with the overall level of radiation - including
    ///    direct and diffuse shading - is separately accounted for).
    ///
    ///    To set a reference point, it was assumed that, in the absence of a reduction due inverter efficiency, output would
    ///    be proportional to the incident solar radiation and therefore that the normalised output would be equal to the
    ///    shading factor - e.g. if 75% of the radiation was transmitted we would get 75% of the power output.
    ///
    ///    The difference between this reference output and the output predicted by the equations representing the actual data
    ///    was assumed to be due to the efficiency reduction of the inverter induced by overshading. The ratio of the 'expected'
    ///    output and the 'actual' output was calculated. This step disaggregates the reduction due to inverter efficiency from
    ///    the reduction simply due to lower incident radiation, to avoid this being double counted. This is tehrefore the
    ///    correction factor needed in HEM to take into consideration the reduced efficiency of inverters when PV panels are
    ///    partially overshaded.
    ///
    ///    2nd order polynomial curves were fitted through the resulting inverter efficiency factors (as a function of the direct
    ///    factor f_sh_dir, resulting in the equations used below). A two stage polynomial is needed for each inverter type
    ///    because of the step change in behaviour when the 'lower limit' referred to in the source is reached.
    pub(crate) fn inverter_efficiency_lookup(&self, f_sh_dir: f64) -> f64 {
        // Calculate inverter efficiency based on direct shading factor and inverter type
        let x = f_sh_dir;
        let inverter_type = self.inverter.r#type();

        match inverter_type {
            InverterType::StringInverter => {
                let thresh = 0.7;
                let a = 2.7666;
                let b = -4.3397;
                let c = 2.2201;
                let d = -1.9012;
                let e = 4.8821;
                let f = -1.9926;
                if x < thresh {
                    1.0f64.min(a * x.powi(2) + b * x + c)
                } else {
                    1.0f64.min(d * x.powi(2) + e * x + f)
                }
            }
            InverterType::OptimisedInverter => {
                let thresh = 0.42;
                let a = 2.7666;
                let b = -4.3397;
                let c = 2.2201;
                let d = -0.2024;
                let e = 0.4284;
                let f = 0.7721;
                if x < thresh {
                    1.0f64.min(a * x.powi(2) + b * x + c)
                } else {
                    1.0f64.min(d * x.powi(2) + e * x + f)
                }
            }
        }
    }
}

impl OnSiteGeneration for PhotovoltaicSystem {
    /// Produce electrical energy (in kWh) from the PV system
    /// according to BS EN 15316-4-3:2017
    #[allow(unused_assignments)] // Python overwrites weighted_f_sh_dir without using it (reported to DESNZ 2025-11-21) - remove this attribute when corrected in future bug fix
    fn produce_energy(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        // Calculate a weighted f_sh_dir from the weighted average of the
        // direct shading factor and unshaded output of each panel
        let mut total_unshaded_energy_produced = 0.;
        let mut weighted_f_sh_dir = 0.;

        for panel in &self.panels {
            let energy_produced = panel.produce_energy(
                self.clone(),
                &self.external_conditions,
                simulation_time_iteration,
                1.0,
            )?;
            total_unshaded_energy_produced += energy_produced;
            weighted_f_sh_dir += panel
                .shading_factors_direct(&self.external_conditions, simulation_time_iteration)?
                * energy_produced;
        }

        if total_unshaded_energy_produced <= 0. {
            weighted_f_sh_dir = 1.;
        } else {
            weighted_f_sh_dir /= total_unshaded_energy_produced;
        }

        let mut total_energy_produced = 0.;

        for panel in &self.panels {
            weighted_f_sh_dir = panel
                .shading_factors_direct(&self.external_conditions, simulation_time_iteration)?;
            let energy_produced = panel.produce_energy(
                self.clone(),
                &self.external_conditions,
                simulation_time_iteration,
                weighted_f_sh_dir,
            )?;
            total_energy_produced += energy_produced;
        }

        // Add energy produced to the applicable energy supply connection via inverter (this will reduce demand)
        // calculate total power going to inverter
        let power_input_inverter = total_energy_produced / simulation_time_iteration.timestep;
        // Inverter will calculate the efficiency and apply ac/dc limits and return energy delivered to the supply connection
        let total_energy_delivered = self
            .inverter
            .produce_energy(power_input_inverter, simulation_time_iteration.index)?;
        let total_energy_lost = total_energy_produced - total_energy_delivered;

        Ok((total_energy_delivered, total_energy_lost))
    }

    /// Return whether this unit is considered inside the building or not
    fn inverter_is_inside(&self) -> bool {
        self.inverter.is_inside()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energy_supply::energy_supply::{
        EnergySupply, EnergySupplyBuilder, EnergySupplyConnection,
    };
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::hem_core::simulation_time::SimulationTimeIterator;
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    fn simulation_time_iterator() -> SimulationTimeIterator {
        simulation_time().iter()
    }

    #[fixture]
    fn external_conditions(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> Arc<ExternalConditions> {
        Arc::new(ExternalConditions::new(
            &simulation_time_iterator,
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
                    start: 180.,
                    end: 135.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: 135.,
                    end: 90.,
                    shading_objects: vec![ShadingObject {
                        object_type: ShadingObjectType::Overhang,
                        height: 2.2,
                        distance: 6.,
                    }],
                },
                ShadingSegment {
                    start: 90.,
                    end: 45.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: 45.,
                    end: 0.,
                    shading_objects: vec![
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
                    ],
                },
                ShadingSegment {
                    start: 0.,
                    end: -45.,
                    shading_objects: vec![ShadingObject {
                        object_type: ShadingObjectType::Obstacle,
                        height: 3.,
                        distance: 8.,
                    }],
                },
                ShadingSegment {
                    start: -45.,
                    end: -90.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: -90.,
                    end: -135.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: -135.,
                    end: -180.,
                    ..Default::default()
                },
            ]
            .into(),
        ))
    }

    fn energy_supply(simulation_time: SimulationTime) -> Arc<RwLock<EnergySupply>> {
        Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ))
    }

    fn energy_supply_connection(simulation_time: SimulationTime) -> EnergySupplyConnection {
        let energy_supply = energy_supply(simulation_time);
        EnergySupply::connection(energy_supply, "pv generation without shading").unwrap()
    }

    fn inverter(
        energy_supply_connection: EnergySupplyConnection,
        inverter_is_inside: bool,
        inverter_peak_power_ac: Option<f64>,
        inverter_type: Option<InverterType>,
        simulation_time: SimulationTime,
    ) -> Inverter {
        Inverter::new(
            energy_supply_connection,
            simulation_time.iter(),
            2.5,
            inverter_peak_power_ac.unwrap_or(0.05),
            inverter_is_inside,
            inverter_type.unwrap_or(InverterType::OptimisedInverter),
        )
    }

    #[fixture]
    fn pv_system_with_energy_supply(
        external_conditions: Arc<ExternalConditions>,
        simulation_time: SimulationTime,
    ) -> (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>) {
        let panels = vec![PhotovoltaicPanel::new(
            2.5,
            PhotovoltaicVentilationStrategy::ModeratelyVentilated,
            30.,
            0.,
            10.,
            2.,
            3.,
            simulation_time.step,
            vec![],
        )];
        let energy_supply = energy_supply(simulation_time);
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "pv generation without shading")
                .unwrap();
        let inverter = inverter(energy_supply_connection, false, None, None, simulation_time);
        let pv = PhotovoltaicSystem::new(external_conditions, panels, inverter);
        (pv, energy_supply)
    }

    #[fixture]
    fn pv_system_with_shading_with_energy_supply(
        external_conditions: Arc<ExternalConditions>,
        simulation_time: SimulationTime,
    ) -> (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>) {
        let panels = vec![PhotovoltaicPanel::new(
            2.5,
            PhotovoltaicVentilationStrategy::ModeratelyVentilated,
            30.,
            0.,
            10.,
            2.,
            3.,
            simulation_time.step,
            vec![
                WindowShadingObject::Overhang {
                    depth: 0.5,
                    distance: 0.5,
                },
                WindowShadingObject::SideFinLeft {
                    depth: 0.25,
                    distance: 0.1,
                },
                WindowShadingObject::SideFinRight {
                    depth: 0.25,
                    distance: 0.1,
                },
            ],
        )];

        let energy_supply = energy_supply(simulation_time);
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "pv generation with shading").unwrap();
        let inverter = inverter(
            energy_supply_connection,
            true,
            Some(0.02),
            Some(InverterType::StringInverter),
            simulation_time,
        );
        let pv = PhotovoltaicSystem::new(external_conditions, panels, inverter);
        (pv, energy_supply)
    }

    #[rstest]
    fn test_is_inside(
        pv_system_with_energy_supply: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        pv_system_with_shading_with_energy_supply: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
    ) {
        let (pv_system, _) = pv_system_with_energy_supply;
        let (pv_system_with_shading, _) = pv_system_with_shading_with_energy_supply;
        assert!(!pv_system.inverter_is_inside());
        assert!(pv_system_with_shading.inverter_is_inside());
    }

    // #[ignore]
    #[rstest]
    fn test_produce_energy(
        pv_system_with_energy_supply: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let (pv_system, energy_supply) = pv_system_with_energy_supply;
        let expected_generation_results = [
            -0.002911179810082315,
            -0.01585915973389526,
            -0.03631681332778666,
            -0.0462218185635626,
            -0.05,
            -0.03841528069730012,
            -0.019985927280177524,
            -0.014819433057321862,
        ];
        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            pv_system.produce_energy(t_it).unwrap();
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation without shading"][t_idx],
                expected_generation_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_produce_energy_multiple_panels(
        external_conditions: Arc<ExternalConditions>,
        simulation_time: SimulationTime,
    ) {
        let energy_supply = energy_supply(simulation_time);
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "pv generation").unwrap();
        let pv_inverter = inverter(
            energy_supply_connection,
            false,
            None,
            Some(InverterType::OptimisedInverter),
            simulation_time,
        );

        let panels = vec![
            PhotovoltaicPanel::new(
                2.5,
                PhotovoltaicVentilationStrategy::ModeratelyVentilated,
                30.,
                0.,
                10.,
                2.,
                3.,
                simulation_time.step,
                vec![],
            ),
            PhotovoltaicPanel::new(
                3.5,
                PhotovoltaicVentilationStrategy::ModeratelyVentilated,
                32.,
                0.,
                12.,
                2.,
                3.,
                simulation_time.step,
                vec![],
            ),
        ];
        let pv_system = PhotovoltaicSystem::new(external_conditions, panels, pv_inverter);

        let expected_results = [
            -0.015826768143293042,
            -0.05,
            -0.05,
            -0.05,
            -0.05,
            -0.05,
            -0.05,
            -0.05,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            pv_system.produce_energy(t_it).unwrap();

            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation"][t_idx],
                expected_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_produce_energy_with_shading(
        pv_system_with_shading_with_energy_supply: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let (pv_system_with_shading, energy_supply) = pv_system_with_shading_with_energy_supply;
        let expected_generation_results = [
            -0.0015507260447403823,
            -0.008732593505451556,
            -0.02,
            -0.02,
            -0.02,
            -0.013971934370722397,
            -0.006779589823217942,
            -0.007020372160065822,
        ];
        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            pv_system_with_shading.produce_energy(t_it).unwrap();
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation with shading"][t_idx],
                expected_generation_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_energy_produced_and_energy_lost(
        pv_system_with_energy_supply: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let (pv_system, _) = pv_system_with_energy_supply;
        let expected_energy_produced = [
            0.002911179810082315,
            0.01585915973389526,
            0.03631681332778666,
            0.0462218185635626,
            0.05,
            0.03841528069730012,
            0.019985927280177524,
            0.014819433057321862,
        ];
        let expected_energy_lost = [
            0.012982394066666526,
            0.022261865114937766,
            0.024010573568482414,
            0.023376882776544372,
            0.02560556391283768,
            0.02392305675244094,
            0.023189444399645136,
            0.021949092776458276,
        ];
        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let (energy_produced, energy_lost) = pv_system.produce_energy(t_it).unwrap();
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

    #[rstest]
    fn test_inverter_efficiency_lookup_for_optimised_inverter(
        pv_system_with_energy_supply: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
    ) {
        let (pv_system, _) = pv_system_with_energy_supply;
        assert_relative_eq!(pv_system.inverter_efficiency_lookup(0.9), 0.993716);
        assert_relative_eq!(pv_system.inverter_efficiency_lookup(0.3), 1.);
    }

    #[rstest]
    fn test_inverter_efficiency_lookup_for_string_inverter(
        external_conditions: Arc<ExternalConditions>,
        simulation_time: SimulationTime,
    ) {
        let pv_string_inverter = inverter(
            energy_supply_connection(simulation_time),
            false,
            None,
            Some(InverterType::StringInverter),
            simulation_time,
        );
        let pv_system_with_string_inverter =
            PhotovoltaicSystem::new(external_conditions, vec![], pv_string_inverter);

        assert_relative_eq!(
            pv_system_with_string_inverter.inverter_efficiency_lookup(0.9),
            0.8613180000000007
        );

        assert_relative_eq!(
            pv_system_with_string_inverter.inverter_efficiency_lookup(0.5),
            0.7419000000000002
        );
    }

    #[rstest]
    fn test_product_energy_zero_ratio_of_rated_output(
        external_conditions: Arc<ExternalConditions>,
        simulation_time: SimulationTime,
    ) {
        let pv_inverter = inverter(
            energy_supply_connection(simulation_time),
            false,
            Some(0.02),
            Some(InverterType::StringInverter),
            simulation_time,
        );

        // Test that ratio_of_rated_output of 0 returns 0 energy
        let panel = PhotovoltaicPanel::new(
            0.,
            PhotovoltaicVentilationStrategy::ModeratelyVentilated,
            30.,
            0.,
            10.,
            2.,
            3.,
            simulation_time.step,
            vec![],
        );
        let pv_system = PhotovoltaicSystem::new(external_conditions, vec![panel], pv_inverter);
        let simulation_time_iteration = simulation_time.iter().next().unwrap();
        assert_eq!(
            pv_system.produce_energy(simulation_time_iteration).unwrap(),
            (0., 0.)
        );
    }
}
