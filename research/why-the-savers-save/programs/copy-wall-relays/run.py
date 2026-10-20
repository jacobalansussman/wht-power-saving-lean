"""run.py (d-relay) -- run relaysat.py cases, each in its own process with a hard wall-clock limit.
usage: python3 -I run.py LIMIT_S 'INST r L MODE [--gauge] [--save f]' ...   (prints one line per case)"""
import sys, subprocess, time
VENV = __import__('sys').executable
HERE = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.')
lim = float(sys.argv[1])
for case in sys.argv[2:]:
    t0 = time.time()
    try:
        p = subprocess.run([VENV, '-I', HERE + '/relaysat.py'] + case.split(), capture_output=True, text=True, timeout=lim)
        print(p.stdout.strip() or p.stderr.strip()[-400:], flush=True)
    except subprocess.TimeoutExpired:
        print('%-60s UNDECIDED (killed at %.0fs)' % (case, time.time() - t0), flush=True)
