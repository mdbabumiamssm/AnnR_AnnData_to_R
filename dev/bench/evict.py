"""Evict a file's pages from the OS page cache (no root needed)."""
import os, sys
for p in sys.argv[1:]:
    fd = os.open(p, os.O_RDONLY)
    os.fsync(fd) if False else None
    os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    os.close(fd)
