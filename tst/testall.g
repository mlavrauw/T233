#
# T233: Algorithms for tensors.
#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage( "T233" );

exclude := [];

TestDirectory(DirectoriesPackageLibrary( "T233", "tst" ),
    rec(
      exitGAP := true,
      exclude := exclude,
      testOptions := rec(compareFunction := "uptowhitespace")
      #rewriteToFile := true,  # enable this line to update tests
    ));
FORCE_QUIT_GAP(1);
