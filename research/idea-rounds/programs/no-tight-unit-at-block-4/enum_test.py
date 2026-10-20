import sys, collections
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b4_admissible.py')).read().split("if __name__ == '__main__':")[0]
exec(src)
reps = enumerate_depth(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
print("frames by depth:", sorted(collections.Counter(depth(f) for f in reps).items()), "expected 1, 30, 280, 960, 1024")
print("K frame depth:", depth(frame(F.K_WORD)))
