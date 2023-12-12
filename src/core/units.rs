pub const JOULES_PER_KILOWATT_HOUR: u32 = 3_600_000;
pub const JOULES_PER_KILOJOULE: u32 = 1000;
pub const WATTS_PER_KILOWATT: u32 = 1000;
pub const LITRES_PER_CUBIC_METRE: u32 = 1000;
pub const SECONDS_PER_HOUR: u32 = 3600;
pub const DAYS_IN_MONTH: [u32; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

pub fn average_monthly_to_annual(list_monthly_averages: [f64; 12]) -> f64 {
    list_monthly_averages
        .iter()
        .enumerate()
        .map(|(month_idx, month_ave)| month_ave * DAYS_IN_MONTH[month_idx] as f64)
        .sum::<f64>()
        / DAYS_IN_MONTH.iter().sum::<u32>() as f64
}

pub fn kelvin_to_celsius(temp_k: f64) -> f64 {
    assert!(temp_k >= 0.0);
    temp_k - 273.15
}

pub fn celsius_to_kelvin(temp_c: f64) -> f64 {
    assert!(temp_c >= -273.15);
    temp_c + 273.15
}
