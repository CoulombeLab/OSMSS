SAINTexpress goes in this folder.

OSMSS uses SAINTexpress 3.6.3 (intensity mode) to score affinity-purification
(AP-MS) data. SAINTexpress is not included with OSMSS. Download it from
https://sourceforge.net/projects/saint-apms/files/ and arrange it like this:

  SAINT/
    SAINTexpress_v3.6.3__2018-03-09/
      bin/
        SAINTexpress-int        (Windows: SAINTexpress-int.exe)

Windows: copy Precompiled_binaries/Windows64/SAINTexpress-int.exe into bin/.
Linux:   copy Precompiled_binaries/Linux64/SAINTexpress-int into bin/ and run
         chmod +x on it, or run "make -j" in the folder to build into bin/.
macOS:   no precompiled binary; run "make -j" in the folder to build into bin/.

Restart OSMSS afterwards. Full instructions: see "SAINTexpress (AP-MS)" in the
main README.md.
