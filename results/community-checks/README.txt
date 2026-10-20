Logs of the community's own checks of the circuit, re-run by the agents on a Linux machine (Colab, Python 3.12.3
with the package versions #3945 pins), 19 October 2026. The community's code was run unchanged, from the pull
requests' own files; nothing here was run on a Mac.

  full3945.out        #3945's complete stated verification, admission/full_replay.py --jobs 3, from scratch with
                      fresh native replays: PASS_TERMINAL_WIDTH5_COMPLETE_CONDITIONAL_RELEASE, 119 of 119 admission
                      controls and 12 of 12 bit controls rejected, 1 h 08 min, 6.2 GB at most
  native-base.out,    #3945's five-stage native replay of both words (native/run_native.py):
  native-extra.out    PASS_COMPLETE_T74_NATIVE_REPLAY, including the 26 widened terminals and the width-5 helper phase
  op-*.out            the exact eight-stream operator replay (ChaloupeMusic's gist, on the operator model of #3187
                      and #3296) over all 360 owners: 0 differences; its controls reject
  reg3945/            #3945's eight registered commands and the community's package runner: all pass
  full3815.out,       #3815 (the fallback): its complete admission with fresh natives, and its registered checks
  reg3815/
  reg4076/, heads-4076.json, lever-4076.txt
                      the draft #4076 that appeared during the build: its registered package and fresh natives pass;
                      the head test counts 57,798 heads on its fill; the half cut would give 1.7504268e-3 (not built)
