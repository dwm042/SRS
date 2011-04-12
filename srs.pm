#!/usr/bin/perl

=head1 NAME

 srs.pm - simple ranking system, coded three different ways.

=head1 SYNOPSIS

 Implements the simple ranking system, as described by Doug Drinen, three different ways.
 Please note that all implementations other than the iterative one tend to have issues
 as many matrices that result from the equations are singular or near singular.

=head1 VERSION

 author   dwmyers
 date     4/06/2011
 modified 4/12/2011

=head1 DESCRIPTION

 srs.pm - simple ranking system, coded three different ways. Note that the simultaneous
 use of Math::MatrixReal and PDL limits the "intelligence" of these choices. This module
 is more for illustration than hard core use.

 These subroutines expect a %teams hash to be passed to them. The keys of the hash should
 be the team names. The required fields of the hash are..

 $teams{NAME}{games_played}
 $teams{NAME}{point_spread}  # sum total point spread for all games played
 $teams{NAME}{played}        # pointer to an array of NAME values, generally loaded like
                             # this:
                             # push @{$teams{NAME}{played}}, NEW_NAME;

 If a matrix is involved, these routines return a 1 or 0. Sorry, guys, but SRS often yields
 a singular or near singular matrix. You then fall back to the iterative solution as described
 by Doug.

 SRS described here:

 http://www.pro-football-reference.com/blog/?p=37

=head1 COPYRIGHT

    This program is copyrighted 2011 by David Myers. All rights are
    reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

package srs;

use warnings;
use strict;
use Math::MatrixReal;
use List::Util;
use PDL;

require Exporter;
our @ISA        = qw(Exporter);
our @EXPORT     = qw(calc_srs calc_srs_lr calc_srs_svd);
our $VERSION    = 0.90;

my $debug = 0;
my $num_teams = 32;

sub calc_srs_svd {
    my $tptr = shift;
    my $correct = shift || 1;
    my $index = 0;
    my %team_map;
    my @mov=();
    for ( keys %$tptr ) {
        $team_map{$_} = $index;
        $tptr->{$_}{mov} = $tptr->{$_}{point_spread}/$tptr->{$_}{games_played};
        $mov[$index] = $tptr->{$_}{mov};
        $index++;
    }
    push @mov, 0;
    my $A = zeroes $num_teams, $num_teams + 1;
    for ( keys %$tptr ) {
        for my $opp ( @{$tptr->{$_}{played}} ) {
            my $i = $team_map{$_};
            my $j = $team_map{$opp};
            $A->set($j, $i, $A->at($j, $i) - 1.0/$tptr->{$_}{games_played});
        }
    }
    for ( 0 .. ( $index - 1 ) ) {
        $A->set( $_, $index, 1.0 );
    }
    my $MOV = pdl @mov;
    my ( $R1, $S, $R2 ) = svd $A;
    if (  abs (prod $S ) < 1.0E-40  ) {
        # singular matrix. Should fall back to iterative solution.
        return 0;
    }
    my $MAT_S = inv stretcher $S;
    my $SRS = $R2->transpose x $MAT_S x $R1->transpose  x $MOV->transpose;
    for ( keys %$tptr ) {
        $tptr->{$_}{srs} = $SRS->at(0, $team_map{$_});
        $tptr->{$_}{sos} = $tptr->{$_}{srs} - $tptr->{$_}{mov} 
    }
    srs_correction( $tptr ) if $correct;
    return 1;

}

sub calc_srs_lr {
    my $tptr = shift;
    my $correct = shift || 1;
    my $index = 1;
    my %team_map;
    my $mov = Math::MatrixReal->new( $num_teams, 1);
    for ( keys %$tptr ) {
        $team_map{$_} = $index;
        $tptr->{$_}{mov} = $tptr->{$_}{point_spread}/$tptr->{$_}{games_played};
        $mov->assign($index, 1, $tptr->{$_}{mov});
        $index++;
    }
    my $matrix = Math::MatrixReal->new( $num_teams, $num_teams );
    $matrix->one;
    for ( keys %$tptr ) {
        for my $opp ( @{$tptr->{$_}{played}} ) {
            my $i = $team_map{$_};
            my $j = $team_map{$opp};
            $matrix->assign($i, $j, $matrix->element($i,$j) - 1.0/$tptr->{$_}{games_played});
        }
    }
    print $matrix if $debug;
    my $lr = $matrix->decompose_LR();
    if ( my ( $dim, $x, $basis ) = $lr->solve_LR( $mov ) ) {
        for ( keys %team_map ) {
            $tptr->{$_}{srs} = $x->element($team_map{$_}, 1 );
            $tptr->{$_}{sos} = $tptr->{$_}{srs} - $tptr->{$_}{mov}
        }
    }
    else {
        # singular matrix. Should fall back to iterative solution.
        return 0;
    }
    srs_correction( $tptr ) if $correct;
    return 1;
}

sub calc_srs {
    my $tptr = shift;
    my $correct = shift || 1;
    for ( keys %$tptr ) {
        $tptr->{$_}{mov} = $tptr->{$_}{point_spread}/$tptr->{$_}{games_played};
        $tptr->{$_}{srs} = $tptr->{$_}{mov};
        $tptr->{$_}{oldsrs} = $tptr->{$_}{srs};
        $tptr->{$_}{sos} = 0;
    }
    my $delta = 10.0;
    my $iters = 0;
    while ( $delta > 0.001 ) {
        $iters++;
        $delta = 0.0;
        for ( keys %$tptr ) {
            my $sos = 0.0;
            for my $g ( @{$tptr->{$_}{played}} ) {
                $sos += $tptr->{$g}{srs};
            }
            $sos /= $tptr->{$_}{games_played};
            $tptr->{$_}{srs} = $tptr->{$_}{mov} + $sos;
            my $newdelt = abs( $sos - $tptr->{$_}{sos} );
            $tptr->{$_}{sos} = $sos;
            $delta = List::Util::max ( $newdelt, $delta );
        }
        for ( keys %$tptr ) {
            $tptr->{$_}{oldsrs} = $tptr->{$_}{srs};
        }
    }
    srs_correction ( $tptr ) if $correct;
    print "iters = $iters\n" if $debug;
    return;
}

# Since the solutions to the srs equation wander all
# over the place 
# any solution SRS = MOV + SOS has an equally valid solution
# SRS + c = MOV + SOS + c
# you have to correct for that by setting the sum 
# of all srs values to average to 0.0.
#
sub srs_correction {
    my $tptr = shift;
    my $sum =  0.0;
    for ( keys %$tptr ) {
        $sum += $tptr->{$_}{srs};
    }
    $sum /= $num_teams;
    for ( keys %$tptr ) {
       $tptr->{$_}{srs} -= $sum;
       $tptr->{$_}{sos} -= $sum;
    }
}
