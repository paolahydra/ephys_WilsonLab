# ephys_WilsonLab
code used to run patching experiments with bilateral piezo stimulation for Okubo et al 2020

runExperiment.m is the main script calling all the specialized workflow to run a patch clamp experiment with calibrated bilateral antennal stimulation.

In particular,
run_continuousAcquisition.m --> recordExperiment.m set the nidaq for continuous input and output streaming and acquisition.

Cameras were also controlled through matlab api (see basic_camera_setup.m). However, they were not actually in use during recording, but only before and after recording (the pc would have crashed. Consider using a separate pc).

