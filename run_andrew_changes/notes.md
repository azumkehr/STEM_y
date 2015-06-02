I got stem to run through to completion for the larger domain!

changes I made to the code and the runfile:

- created directory /home/ecampbell_lab/COS/fog/output/. It didn't
  exist; this caused STEM to crash because STEM was told to put output
  there (in the run file)

- changed starting time (isthr) from 0 to 1 and starting date from
  2015 05 04 (day 124) to 2015 05 05 (day 125) in the run file.
  METEO3D starts at 01:00:00 on 2015 day 125.  STEM crashed when the
  run file told it to start at 00:00:00 on day 124.

- changed duration of run (iperiod) from 12000-some hours to 24.
  METEO3D only contains 24 hours of data.

- READ_IOAPI (and, by inheritance, my spinoff
  READ_IOAPI_FIND_PREVIOUS) are written to fail and exit if their
  target file does not exist.  I changed this behavior to issue a
  warning and continue.

- increased mxgr from 400 to 800 in two subroutines within
  output_ioapi.f: subroutine prtemp and subroutine prtout.  This is
  needed to accomodate the larger domain.  mxgr seems to be used to
  calculate the maximum number of records that may be written to an
  output file.  If it tries to write a larger output array than that
  maximum it triggers a fortran STOP.  I'm not sure of two things: (1)
  why this maximum is in there -- maybe to prevent a runaway writeout
  if some erroneous STEM setup tries to produce a huge output file?
  And (2) why it's hard-coded instead of calculated from the domain
  size or something.  Elliott, do you know anything about this?  I'm
  tempted to change it along the lines of (2).  Right now I've just
  increased the hard-coded value from 400 to 800.

And presto, STEM ran through!

-Tim
