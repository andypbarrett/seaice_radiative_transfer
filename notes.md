Portland Group Compiler Installation
https://www.pgroup.com/resources/docs/17.10/x86/pgi-install-guide/index.htm#install-linux-pgi

## From Bonnie on tests
hi Andy,
This is great. Yep, the direct/diffuse weightings looks good.
As for strategies going forward-- I agree best to check numbers-- bummer
about formats (when are they ever the same though??!). I'd expect checks of
albedo, transmittance, and F_sfc_vs and F_sfc_ir to be a comprehensive check.
This is great progress!
Bonnie

On Mon, Feb 28, 2022 at 6:03 PM Andy Barrett <Andrew.Barrett@colorado.edu> wrote:

    Hi Bonnie

    Good to see you today.  I'm attaching a plot showing the same data as in
    Figure 10 of Breigleib and Light - model spectral albedo for bare ice
    and ponds.  It doesn't show the observed spectral albedos because I
    don't have that data.

    I calculated weighted averages of direct and diffuse albedos using the
    fractions of direct and diffuse visible and nir albedos, using nir
    direct and diffuse fractions for model bands 2 and 3, and vis direct and
    diffuse fractions for band 1.  Is this correct?

    I've attached sample output to show these variables.  The vs and nir
    fractions are in the surface absorption and albedos sections (~ Line 153).

    I'm thinking about testing strategies going forward.  To check the code,
    it makes more sense to me to check numbers.  I would just do a diff for
    the .dat files in the zipfile you sent ad the files I produce but the
    formats are slightly different, although the numbers look similar.  So
    my question is, are there some key outputs to check (albedos,
    transmittance) or should I do the whole lot.


