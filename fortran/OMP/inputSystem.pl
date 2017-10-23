#!/usr/bin/perl

# enter kinetic law
print "Kinetic law: "; #" (combinatorial(c) or mass action (m)): ";
$kinetic_law = <STDIN>;
chomp($kinetic_law);


# enter reactions
@reactions = ();
$R_ind = 0;
@reversible_ind = ();
print "\nReactions\n";
print "R[0]: ";
while ( <STDIN> )
{
    # remove \n from end of $_
    chomp;
    if (/\W$/)
    {
        ( $_ =~ /\.$/ ) ? ( $stop = 1 ) : " ";
        $_ =~ s/\W$//;
    }

    if ( /<->/ )
    {
        push @reversible_ind, $R_ind;
        @reactants_products = split(/<->/, $_);
        push @reactions, "@reactants_products[0]"."->"."@reactants_products[1]";
        push @reactions, "@reactants_products[1]"."->"."@reactants_products[0]";
        $R_ind+=2;
    }
    else
    {
        push @reactions, $_;
        $R_ind+=1;
    }

    ($stop) ? (last) : print "R[$R_ind]: ";
}
print "\n\n";
$num_reactions = @reactions;


# identify species
$species = "";
for (@reactions)
{
    @reactants_products = split(/->/, $_);
    @reactants = split(/\+/, @reactants_products[0]);
    @products  = split(/\+/, @reactants_products[1]);
    
    # trim whitespaces
    for (@reactants)
    {
        $_ =~ s/^\s*(.*?)\s*$/$1/;
    }
    for (@products)
    {
        $_ =~ s/^\s*(.*?)\s*$/$1/;
    }

    # find all new species
    for (@reactants)
    {
        my @coef_species = split(/\s/,$_);
        my $len = split(/\s/,$_);

        ( $len == 2 ) ? ( ( $species !~ /@coef_species[1]/ ) ? ( $species = $species.@coef_species[1]." " ) : " " )
                      : ( ( $species !~ /@coef_species[0]/ ) ? ( $species = $species.@coef_species[0]." " ) : " " );
    }

    for (@products)
    {
        my @coef_species = split(/\s/,$_);
        my $len = split(/\s/,$_);

        ( $len == 2 ) ? ( ( $species !~ /@coef_species[1]/ ) ? ( $species = $species.@coef_species[1]." " ) : " " )
                      : ( ( $species !~ /@coef_species[0]/ ) ? ( $species = $species.@coef_species[0]." " ) : " " );
    }
}
@species = split(/\s/,$species);
$num_species = @species;


# reaction rates
print "Reaction rates:    ";
$_ = <STDIN>;
chomp();
@rates = split(/\s/,$_);
#print "\n";

# reaction rates
print "Initial condition: ";
$_ = <STDIN>;
chomp();
@X0 = split(/\s/,$_);
print "\n";


# build stoichiometry
$nu_in  = "";
$nu_out = "";
$nu     = "";
@prop   = ();
$R_ind  = 0;
for (@reactions)
{
    @reactants_products = split(/->/, $_);
    @reactants = split(/\+/, @reactants_products[0]);
    @products  = split(/\+/, @reactants_products[1]);

    # trim whitespaces
    for (@reactants)
    {
        $_ =~ s/^\s*(.*?)\s*$/$1/;
    }
    for (@products)
    {
        $_ =~ s/^\s*(.*?)\s*$/$1/;
    }


    # initialize
    my @coef_in  = ();
    my @coef_out = ();
    my @coef     = ();
    for $i (0..$num_species-1)
    {
        @coef_in[$i]  = 0;
        @coef_out[$i] = 0;
    }

    for (@reactants)
    {
        my @coef_species = split(/\s/,$_);
        my $len = split(/\s/,$_);


        if ( $len == 2 )
        {
            $c = @coef_species[0];
            $s = @coef_species[1];
        }
        else
        {
            $c = 1;
            $s = @coef_species[0];
        }

        $ind = 0;
        for (@species)
        {
            if ( /$s/ )
            {
                @coef_in[$ind] = $c;
            } 
            $ind++;
        }
    }

    for (@products)
    {
        my @coef_species = split(/\s/,$_);
        my $len = split(/\s/,$_);


        if ( $len == 2 )
        {
            $c = @coef_species[0];
            $s = @coef_species[1];
        }
        else
        {
            $c = 1;
            $s = @coef_species[0];
        }

        $ind = 0;
        for (@species)
        {
            if ( /$s/ )
            {
                @coef_out[$ind] = $c;
            } 
            $ind++;
        }
    }

    for $i (0..$num_species-1)
    {
        @coef[$i] = @coef_out[$i] - @coef_in[$i];
    }
    push @nu_in,  @coef_in;
    push @nu_out, @coef_out;
    push @nu,     @coef;
    

    # propensities
    $tmp = $R_ind + 1;
    @prop[$R_ind] = "c($tmp)";
    if ( $kinetic_law eq "combinatorial" )
    {
        for $s (1..$num_species)
        {
            if ( @coef_in[$s-1] >= 1 )
            {
                @prop[$R_ind] = @prop[$R_ind]." * X($s)";
            }
            if ( @coef_in[$s-1] > 1 )
            {
                for $j (1..@coef_in[$s-1]-1)
                {
                    @prop[$R_ind] = @prop[$R_ind]." * ( X($s) - $j )";
                }
            }
        }
    }
    if ( $kinetic_law eq "mass action" )
    {
        for $s (1..$num_species)
        {
            if ( @coef_in[$s-1] == 1 )
            {
                @prop[$R_ind] = @prop[$R_ind]." * X($s)";
            }
            if ( @coef_in[$s-1] > 1 )
            {
                @prop[$R_ind] = @prop[$R_ind]." * X($s)**@coef_in[$s-1]";
            }
        }
    }


    # Jacobian
    @propJac[$R_ind] = [()];
    if ( $kinetic_law eq "combinatorial" )
    {
        for $s (1..$num_species)
        {
            $propJac[$R_ind][$s-1] = "c($tmp)";

            for $s1 (1..$num_species)
            {
                if ($s1 != $s)
                {
                    if ( @coef_in[$s1-1] >= 1 )
                    {
                        # print "$s1\n";
                        # print "@coef_in[$s1-1]\n";
                        $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * X($s1)";
                    }
                    if ( @coef_in[$s1-1] > 1 )
                    {
                        # print "@coef_in[$s1-1]\n";
                        for $j (1..@coef_in[$s1-1]-1)
                        {
                            $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * ( X($s1) - $j )";
                        }
                    }
                }
            }

            ( @coef_in[$s-1] == 0 ) ? ( $propJac[$R_ind][$s-1] = "0" ) : " " ;
            ( @coef_in[$s-1] == 1 ) ? " " : " " ;
            ( @coef_in[$s-1] == 2 ) ? ( $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * ( 2*X($s) - 1 )" ) : " " ;
            ( @coef_in[$s-1] == 3 ) ? ( $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * ( 3*X($s)**2 - 6*X($s) + 2)" ) : " " ;
            ( @coef_in[$s-1] == 4 ) ? ( $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * ( 4*X($s)**3 - 18*X($s)**2 + 22*X($s) -6)" ) : " " ;
        }
    }
    if ( $kinetic_law eq "mass action" )
    {
        for $s (1..$num_species)
        {
            $propJac[$R_ind][$s-1] = "c($tmp)";

            for $s1 (1..$num_species)
            {
                if ($s1 != $s)
                {
                    if ( @coef_in[$s1-1] == 1 )
                    {
                        $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * X($s1)";
                    }
                    if ( @coef_in[$s1-1] > 1 )
                    {
                        $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * X($s1)**@coef_in[$s1-1]";
                    }
                }
            }

            ( @coef_in[$s-1] == 0 ) ? ( $propJac[$R_ind][$s-1] = "0" ) : " " ;
            ( @coef_in[$s-1] == 1 ) ? " " : " " ;
            ( @coef_in[$s-1] == 2 ) ? ( $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * 2*X($s)" ) : " " ;
            ( @coef_in[$s-1]  > 2 ) ? ( $propJac[$R_ind][$s-1] = $propJac[$R_ind][$s-1]." * @coef_in[$s-1]*X($s)" ) : " " ;
        }
    }

    $R_ind++;
}


# relaxation rates
for $r1 (0..$num_reactions-1)
{
    for $s (0..$num_species-1)  
    {
        $ind = $r1*$num_species + $s;
        if ( $propJac[$r1][$s] ne "0" )
        {
            my $tmp = abs(@nu[$ind]);
            ( abs(@nu[$ind]) == 0 ) ? " " : " ";
            ( abs(@nu[$ind]) == 1 ) ? ( (@nu[$ind]>0) ? ( $relRates[$r1] = $relRates[$r1]." - $propJac[$r1][$s]" ) : ( $relRates[$r1] = $relRates[$r1]." + $propJac[$r1][$s]" ) ) : " ";
            ( abs(@nu[$ind])  > 1 ) ? ( (@nu[$ind]>0) ? ( $relRates[$r1] = $relRates[$r1]." - @nu[$ind]*$propJac[$r1][$s]" ) : ( $relRates[$r1] = $relRates[$r1]." + $tmp*$propJac[$r1][$s]" ) )   : " ";
        }
    }
}
for $r1 (@reversible_ind)
{
    $r2 = $r1+1;
    for $s (0..$num_species-1)  
    {
        $ind = $r2*$num_species + $s;
        if ( $propJac[$r2][$s] ne "0" )
        {
            # print "@nu[$ind]\n";
            my $tmp = abs(@nu[$ind]);
            ( abs(@nu[$ind]) == 0 ) ? " " : " ";
            ( abs(@nu[$ind]) == 1 ) ? ( (@nu[$ind]>0) ? ( $relRates[$r1] = $relRates[$r1]." - $propJac[$r2][$s]" ) : ( $relRates[$r1] = $relRates[$r1]." + $propJac[$r2][$s]" ) )  : " ";
            ( abs(@nu[$ind])  > 1 ) ? ( (@nu[$ind]>0) ? ( $relRates[$r1] = $relRates[$r1]." - @nu[$ind]*$propJac[$r2][$s]" ) : ( $relRates[$r1] = $relRates[$r1]." + $tmp*$propJac[$r2][$s]" ) ) : " ";
        }
        $relRates[$r2] = @relRates[$r1];
    }
}


# 
# @nu_times_X = ();
# for $s (0..$num_species-1)
# {
#     for $r (0..$num_reactions-1)
#     {
#         $ind = $r*$num_species + $s;
#         if ( @nu[$ind] ne "0" )
#         {
#             my $tmp = $r + 1;
#             ( @nu[$ind] eq "1" ) ? : ( @nu_times_X[$s] = " + ".@nu_times_X[$s]."X($tmp)" );
#         }
#     }
# }



open(fh, ">module_ChemicalSystem.f90") or die "Can't open data";

$var = <<"EOF";
module ChemicalSystem
    implicit none

    integer(kind=4), parameter                  :: num_reactions = $num_reactions
    integer(kind=4), parameter                  :: num_species   = $num_species
EOF
print fh "$var";
print fh "    integer(kind=4), parameter, dimension($num_species,$num_reactions)  :: nu = reshape((/";
print fh " @nu[0]";
for $r (1..$num_reactions*$num_species-1)
{
    print fh ", @nu[$r]";
}
print fh " /),(/$num_species,$num_reactions/))\n";

print fh "       real(kind=8), dimension($num_species)               :: X0 = (/ ";
print fh "@X0[0]";
for $i (1..$num_species-1)
{
    print fh ", @X0[$i]";
}
print fh " /)\n";
print fh "       real(kind=8), dimension($num_reactions)               :: c = (/ ";
print fh "@rates[0]";
for $i (1..$num_reactions-1)
{
    print fh ", @rates[$i]";
}
print fh " /)\n";

$var = <<"EOF";

    contains

    function propensities(X) result(a)
        real(kind=8), dimension(num_species), intent(in) :: X
        real(kind=8), dimension(num_reactions)           :: a

EOF
print fh "$var";

for $r (1..$num_reactions)
{
    print fh "        a($r) = $prop[$r-1]\n";
}
print fh "    end function\n";

$var = <<"EOF";

    function propJacobian(X) result(Ja)
        real(kind=8), dimension(num_species), intent(in)   :: X
        real(kind=8), dimension(num_reactions,num_species) :: Ja

EOF
print fh "$var";

for $s (1..$num_species)
{
    for $r (1..$num_reactions)
    {
        $r1 = $r - 1;
        $s1 = $s - 1;
        print fh "        Ja($r,$s) = $propJac[$r1][$s1]\n";
    }
}
print fh "    end function\n";



$var = <<"EOF";

    function relaxRates(X) result(lambda)
        real(kind=8), dimension(num_species), intent(in) :: X
        real(kind=8), dimension(num_reactions)           :: lambda

EOF
print fh "$var";

for $r (1..$num_reactions)
{
    print fh "        lambda($r) = $relRates[$r-1]\n";
}
print fh "    end function\n";

print fh "end module ChemicalSystem\n";

close(fh) or die "Can't close file";
