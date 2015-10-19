 Fix outstanding problem when reading in files that have a single time
coordinate. These give cubes which do not concatenate properly with others.

This happens a lot when reading HadGEM2 (or other Met Office) data. This is
because the files start and finish in November. When constraining on time and
requiring complete years this yields cubes with a single time point which do
not concatenate with following cubes.
