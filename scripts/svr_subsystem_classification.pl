# -*- perl -*-
# This is a SAS Component
########################################################################
# Copyright (c) 2003-2008 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
########################################################################

use strict;
use SeedEnv;

=head1 svr_subsystem_classification

Extend a set of subsystems names to classifications

Example:

    svr_subsystem_classification < table.with.ss-name.as.last.column > extended.table

=head2 Output

A table with one added columnns (comma-delimited list of category and sub-category).
Lines in the incoming table that do not match are written to STDERR.

=cut

my $usage = "svr_roles_to_subsys svr_subsystem_classification < table.with.ss-name.as.last.column > extended.table 2> nonmatching.rows";

use Getopt::Long;
my $url    = '';
my $column = undef;
my $rc = GetOptions( "c=i" => \$column,
                     "url=s" => \$url
                   );

$rc or print STDERR $usage and exit;


# Get the server object.
my $sapServer = SAPserver->new(url => $url);

# The main loop processes chunks of input, 1000 lines at a time.
while (my @tuples = ScriptThing::GetBatch(\*STDIN, undef, $column)) {
    # Ask the server for the subsystem classifications
    my $classHash = $sapServer->classification_of(-ids => [ map { $_->[0] } @tuples]);
    
    # Output the results for these roles.
    for my $tuple (@tuples) {
        # Get this line and the subsystem name.
        my ($ss_name, $line) = @$tuple;
	
	# Output the line with the subsystem and classification appended.
	print STDOUT (join("\t", ($line, join("\t", @ { $classHash->{$ss_name} }) )). "\n");
    }
}
exit(0);
