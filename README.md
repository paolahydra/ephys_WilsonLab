# ephys_WilsonLab
code written by Paola Patella and used to run patch clamp experiments combined with bilateral piezo stimulation for Okubo, Patella, D'Alessandro and Wilson (2020).

# custom Matlab experiment controller for setting up and recording patch clamp + sensory stimulation experiments  
runExperiment.m is the main script calling all the specialized workflow to run a patch clamp experiment with calibrated bilateral antennal stimulation.

In particular,
run_continuousAcquisition.m --> recordExperiment.m sets the nidaq for continuous input and output streaming and acquisition.

Cameras were also controlled through matlab api (see basic_camera_setup.m). However, they were not actually in use during recording, but only before and after recording (the pc would have crashed. Consider using a separate pc).

