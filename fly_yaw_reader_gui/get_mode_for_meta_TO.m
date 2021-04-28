function [ mode ] = get_mode_for_meta_TO( mode_meta )
% Convention is mode == the value that ScaledOut represents
% which is quite wrong/counterintuitive.... but hey... that's houw it was
% used here



if mode_meta > 0 && mode_meta < 3.5
    mode = 'V'; % current clamp
elseif mode_meta >= 3.5 && mode_meta < 7
    mode = 'I'; % voltage clamp
elseif mode_meta <= 0
    error('Turn on the amplifier.')
else
    display(['MODE: ',num2str(mode_meta)])
    error('MODE telegraph out of range.')
end

