use crate::compare_floats::{max_of_2, min_of_2};
use crate::input::ElectricBattery as ElectricBatteryInput;
use atomic_float::AtomicF64;
use std::sync::atomic::Ordering;

#[derive(Debug)]
pub struct ElectricBattery {
    /// the maximum capacity of the battery
    capacity: f64,
    /// charge/discharge round trip efficiency of battery system (between 0 & 1)
    charge_discharge_efficiency: f64,
    /// the current energy stored in the battery at the
    /// end of hour (kWh)
    current_energy_stored: AtomicF64,
}

/// Arguments:
/// * `elec_demand` - the supply (-ve) or demand (+ve) to/on the electric battery (kWh)
impl ElectricBattery {
    pub fn new(capacity: f64, charge_discharge_efficiency: f64) -> Self {
        Self {
            capacity,
            charge_discharge_efficiency,
            current_energy_stored: Default::default(),
        }
    }

    pub fn from_input(input: ElectricBatteryInput) -> Self {
        Self::new(input.capacity, input.charge_discharge_efficiency)
    }

    pub fn charge_discharge_battery(&self, elec_demand: f64) -> f64 {
        let energy_available_to_charge_battery = if elec_demand < 0. {
            // Charging battery
            (-elec_demand) * self.charge_discharge_efficiency.powf(0.5)
        } else {
            // Discharging battery
            (-elec_demand) / self.charge_discharge_efficiency.powf(0.5)
        };
        // Charge/discharge the battery by the amount available
        let current_energy_stored_unconstrained =
            self.current_energy_stored.load(Ordering::SeqCst) + energy_available_to_charge_battery;
        let prev_energy_stored = self.current_energy_stored.load(Ordering::SeqCst);
        // Energy stored cannot be > battery capacity or < 0
        self.current_energy_stored.store(
            min_of_2(
                self.capacity,
                max_of_2(0., current_energy_stored_unconstrained),
            ),
            Ordering::SeqCst,
        );
        let energy_accepted_by_battery =
            self.current_energy_stored.load(Ordering::SeqCst) - prev_energy_stored;

        // Return the supply/demand energy the battery can accept (including charging/discharging losses)
        if elec_demand < 0. {
            -energy_accepted_by_battery / self.charge_discharge_efficiency.powf(0.5)
        } else {
            -energy_accepted_by_battery * self.charge_discharge_efficiency.powf(0.5)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_ulps_eq};
    use rstest::*;

    #[fixture]
    pub fn elec_battery() -> ElectricBattery {
        ElectricBattery::new(2000., 0.8)
    }

    #[rstest]
    pub fn test_charge_discharge_battery(elec_battery: ElectricBattery) {
        // supply to battery exceeds limit
        assert_relative_eq!(
            elec_battery.charge_discharge_battery(-1_000_000.),
            -2000. / 0.8f64.powf(0.5),
            max_relative = 1e-7
        );
        // demand on battery exceeds limit
        assert_relative_eq!(
            elec_battery.charge_discharge_battery(100_000_000.),
            2000. * 0.8f64.powf(0.5),
            max_relative = 1e-7
        );
        // normal charge
        assert_ulps_eq!(elec_battery.charge_discharge_battery(-200.), -200.);
        // normal discharge
        assert_ulps_eq!(elec_battery.charge_discharge_battery(100.), 100.);
    }
}
