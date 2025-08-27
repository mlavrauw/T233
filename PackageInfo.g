#
# T233: Algorithms for tensors.
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "T233",
Subtitle := "Algorithms for tensors",
Version := "0.1",
Date := "28/08/2019", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Nour",
    LastName := "Alnajjarine",
    WWWHome := "https://sites.google.com/view/nouralnajjarine/home",
    Email := "nour.alnajjarine@math.uniri.hr",
    PostalAddress := "Nour Alnajjarine, Faculty of Mathematics, University of Rijeka, Croatia",
    Place := "Rijeka, Croatia",
    Institution := "University of Rijeka",
  ),
  rec(
    IsAuthor := true,
    IsMaintainer := true,
    FirstNames := "Michel",
    LastName := "Lavrauw",
    WWWHome := "https://mlavrauw.github.io",
    Email := "michel.lavrauw@upr.si",
    PostalAddress := "Michel Lavrauw, Department of Mathematics, University of Primorska, Koper, Slovenia",
    Place := "Koper, Slovenia",
    Institution := "University of Primorska",
  ),
],







SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/mlavrauw/T233",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := "https://mlavrauw.github.io/T233/",
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=  "",

PackageDoc := rec(
  BookName  := "T233",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Algorithms for tensors",
),

Dependencies := rec(
  GAP := ">= 4.9",
  NeededOtherPackages := ["FinInG", ">=1.5.6"],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

#Keywords := [ "TODO" ],

));
