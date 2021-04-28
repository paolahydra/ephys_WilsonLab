function [currentA, voltageA, currentB, voltageB] = get_dual_scaled_voltage_and_current( trial_data )

[ currentA, voltageA ] = get_scaled_voltage_and_current_A( trial_data );
[ currentB, voltageB ] = get_scaled_voltage_and_current_B( trial_data );

end

