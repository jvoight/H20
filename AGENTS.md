# Working on H20

This repository is the canonical home for the paper's computational code
and associated data. Make subsequent code changes here.

After completing a requested code change, run the relevant checks described
in the directory README, commit the task's changes, and push to the current
branch's upstream. The user has requested this ongoing commit-and-push
workflow. Check for concurrent edits before staging files.

Keep source-data provenance, saved results, and reproduction instructions
together. Magma may return exit status zero after an error; inspect its
output and require the expected PASS messages.

Keep temporary build files, executables, and one-off manuscript editing
scripts out of the repository. The manuscript itself is maintained
separately; update its code references when relevant paths or results change.
