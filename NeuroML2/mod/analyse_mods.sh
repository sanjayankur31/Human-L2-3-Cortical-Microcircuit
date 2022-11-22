# CaDynamics.mod  Ca_LVA.mod  Gfluct.mod  Im.mod   K_T.mod    Nap.mod   NMDA.mod          ProbUDFsyn.mod  tonic.mod
# Ca_HVA.mod      epsp.mod    Ih.mod      K_P.mod  Kv3_1.mod  NaTg.mod  ProbAMPANMDA.mod  SK.mod
#
echo "Compiling mod files"
nrnivmodl .

parallel_run() {
    parallel -j $1 --\
        "pynml-modchananalysis SK -stepV 5 -temperature 6.3 -v"\
        "pynml-modchananalysis Kv3_1 -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis NaTg -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis Nap -stepV 5 -temperature 6.3 -duration 20000"\
        "pynml-modchananalysis K_P -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis K_T -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis Ca_HVA -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis Ca_LVA -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis Ih -stepV 5 -temperature 6.3"\
        "pynml-modchananalysis Im -stepV 5 -temperature 6.3"

}

serial_run () {
    pynml-modchananalysis SK -stepV 5 -temperature 6.3 -v
    pynml-modchananalysis Kv3_1 -stepV 5 -temperature 6.3
    pynml-modchananalysis NaTg -stepV 5 -temperature 6.3
    pynml-modchananalysis Nap -stepV 5 -temperature 6.3 -duration 20000
    pynml-modchananalysis K_P -stepV 5 -temperature 6.3
    pynml-modchananalysis K_T -stepV 5 -temperature 6.3
    pynml-modchananalysis Ca_HVA -stepV 5 -temperature 6.3
    pynml-modchananalysis Ca_LVA -stepV 5 -temperature 6.3
    pynml-modchananalysis Ih -stepV 5 -temperature 6.3
    pynml-modchananalysis Im -stepV 5 -temperature 6.3
}

echo "Analysing mod files"

serial_run

## if gnu parallel is available, use the parallel run, argument is number of parallel jobs:
# parallel_run 6
