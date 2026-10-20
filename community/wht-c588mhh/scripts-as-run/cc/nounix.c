/* nounix (agent c-comparator): run a command with socket(AF_UNIX, ...) forbidden.
   Stand-in for systemd-run --property=RestrictAddressFamilies=~AF_UNIX, which needs systemd.
   Like systemd, it filters the socket() system call only. The filter is inherited by all children. */
#define _GNU_SOURCE
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/prctl.h>
#include <sys/socket.h>
#include <sys/syscall.h>
#include <linux/audit.h>
#include <linux/filter.h>
#include <linux/seccomp.h>

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: nounix cmd [args...]\n"); return 2; }
  struct sock_filter filter[] = {
    /* 0 */ BPF_STMT(BPF_LD | BPF_W | BPF_ABS, offsetof(struct seccomp_data, arch)),
    /* 1 */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, AUDIT_ARCH_X86_64, 1, 0),
    /* 2 */ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_KILL_PROCESS),   /* other architectures (32-bit calls): kill */
    /* 3 */ BPF_STMT(BPF_LD | BPF_W | BPF_ABS, offsetof(struct seccomp_data, nr)),
    /* 4 */ BPF_JUMP(BPF_JMP | BPF_JGE | BPF_K, 0x40000000, 0, 1), /* x32 calls: kill */
    /* 5 */ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_KILL_PROCESS),
    /* 6 */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, __NR_socket, 0, 3),
    /* 7 */ BPF_STMT(BPF_LD | BPF_W | BPF_ABS, offsetof(struct seccomp_data, args[0])),
    /* 8 */ BPF_JUMP(BPF_JMP | BPF_JEQ | BPF_K, AF_UNIX, 0, 1),
    /* 9 */ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ERRNO | (EAFNOSUPPORT & SECCOMP_RET_DATA)),
    /* 10 */ BPF_STMT(BPF_RET | BPF_K, SECCOMP_RET_ALLOW),
  };
  struct sock_fprog prog = { .len = (unsigned short)(sizeof(filter) / sizeof(filter[0])), .filter = filter };
  if (prctl(PR_SET_NO_NEW_PRIVS, 1, 0, 0, 0) != 0) { perror("nounix: PR_SET_NO_NEW_PRIVS"); return 111; }
  if (prctl(PR_SET_SECCOMP, SECCOMP_MODE_FILTER, &prog) != 0) { perror("nounix: PR_SET_SECCOMP"); return 111; }
  execvp(argv[1], argv + 1);
  perror("nounix: execvp");
  return 127;
}
