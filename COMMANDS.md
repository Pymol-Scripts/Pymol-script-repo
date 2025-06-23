

## Draw_Protein_Dimensions.py

### draw_Protein_Dimensions
```None```

### draw_BB
```None```
 

## aKMT_Lys_pred.py

### get_resis_from_resn
```None```

### match_peptides
```None```
 

## aaindex.py

### aaindex2b
```
    DESCRIPTION

        "aaindex" looks up the Amino Acid Index from
          http://www.genome.jp/aaindex/
        for the given key and assignes b-factors to the given selection. Unknown
        residues get the average index value assigned.

    USAGE

        aaindex2b [key [, selection]]

    ARGUMENTS

        key = string: Key of AAindex entry

        selection = string: atoms to assign b-factors {default: (all)}

    EXAMPLE

        # Hydropathy index by Kyte-Doolittle
        aaindex2b KYTJ820101
        spectrumany b, white yellow forest
        show surface
```

### pmf
```
    DESCRIPTION

        Potential of Mean Force

    ARGUMENTS

        key = string: aaindex key

        cutoff = float: distance cutoff {default: 7.0}
        cutoff = (float, float): distance shell

        selection1 = string: atom selection {default: (name CB)}

        selection2 = string: atom selection {default: selection1}

    NOTES

        Does also support a list of keys and a list of cutoffs to deal with
        multiple distance shells.

    EXAMPLES

        # distance dependent c-beta contact potentials
        pmf SIMK990101, 5,         /2x19//A//CB
        pmf SIMK990102, [5, 7.5],  /2x19//A//CB
        pmf [SIMK990101, SIMK990102, SIMK990103], [0, 5, 7.5, 10], /2x19//A//CB

        # interface potential
        sidechaincenters 2x19_scc, 2x19
        pmf KESO980102, 7.0, /2x19_scc//A, /2x19_scc//B
        distance /2x19_scc//A, /2x19_scc//B, cutoff=7.0
```
 

## anglebetweenhelices.py

### helix_orientation
```
    DESCRIPTION

        Get the center and direction of a helix as vectors. Will only work
        for helices and gives slightly different results than loop_orientation.
        Averages direction of C(i)->O(i) bonds.

    USAGE

        helix_orientation selection [, visualize [, sigma_cutoff]]

    ARGUMENTS

        selection = string: atom selection of helix

        visualize = 0 or 1: show fitted vector as arrow {default: 1}

        sigma_cutoff = float: drop outliers outside
        (standard_deviation * sigma_cutoff) {default: 1.5}

    SEE ALSO

        angle_between_helices, helix_orientation_hbond, loop_orientation, cafit_orientation
```

### helix_orientation_hbond
```
    DESCRIPTION

        Get the center and direction of a helix as vectors. Will only work
        for alpha helices and gives slightly different results than
        helix_orientation. Averages direction of O(i)->N(i+4) hydrogen bonds.

    USAGE

        helix_orientation selection [, visualize [, cutoff]]

    ARGUMENTS

        cutoff = float: maximal hydrogen bond distance {default: 3.5}

    SEE ALSO

        helix_orientation
```

### loop_orientation
```
    DESCRIPTION

        Get the center and approximate direction of a peptide. Works for any
        secondary structure.
        Averages direction of N(i)->C(i) pseudo bonds.

    USAGE

        loop_orientation selection [, visualize]

    SEE ALSO

        helix_orientation
```

### cafit_orientation
```
    DESCRIPTION

        Get the center and direction of a peptide by least squares
        linear fit on CA atoms.

    USAGE

        cafit_orientation selection [, visualize]

    NOTES

        Requires python module "numpy".

    SEE ALSO

        helix_orientation
```

### angle_between_helices
```
    DESCRIPTION

        Calculates the angle between two helices

    USAGE

        angle_between_helices selection1, selection2 [, method [, visualize]]

    ARGUMENTS

        selection1 = string: atom selection of first helix

        selection2 = string: atom selection of second helix

        method = string: function to calculate orientation {default: helix_orientation}
                 or int: 0: helix_orientation, 1: helix_orientation_hbond,
                         2: loop_orientation, 3: cafit_orientation

        visualize = 0 or 1: show fitted vector as arrow {default: 1}

    SEE ALSO

        helix_orientation, helix_orientation_hbond, loop_orientation, cafit_orientation
```
 

## annotate_v.py

### annotate_v
```None```
 

## b2transparency.py

### b2transparency
```
    DESCRIPTION

        Set surface (or other) transparency for each atom scaled by b-factor.

        Does not work for all, but for some transparency settings (for example
        transparency, sphere_transparency)

    ARGUMENTS

        selection = string: atom selection {default: all}

        setting = string: setting name {default: transparency}

        minimum = float: b-factor range minimum {default: automatic}

        maximum = float: b-factor range maximum {default: automatic}

        var = string: numeric atomic property like b or q {default: b}

    SEE ALSO

        spectrum, cartoon putty
```
 

## bbPlane.py

### bbPlane
```
    DESCRIPTION

        Draws a plane across the backbone for a selection

    ARGUMENTS

        selection = string: protein object or selection {default: (all)}

        color = string: color name or number {default: white}

        transp = float: transparency component (0.0--1.0) {default: 0.0}

        state = integer: object state, 0 for all states {default: 1}

    NOTES

        You need to pass in an object or selection with at least two
        amino acids.  The plane spans CA_i, O_i, N-H_(i+1), and CA_(i+1)
```
 

## cb_colors.py
 

## ccp4_contact.py

### ccp4_contact
```None```
 

## ccp4_ncont.py

### ccp4_ncont
```None```
 

## ccp4_pisa.py

### ccp4_pisa
```None```
 

## center_of_mass.py

### com
```None```

### get_com
```
    DESCRIPTION

       Calculates the center of mass

       Author: Sean Law
       Michigan State University
       slaw (at) msu . edu
```
 

## centroid.py

### centroid
```None```
 

## cgoCircle.py

### cgoCircle
```
    Create a CGO circle

    PARAMS
          x, y, z
            X, Y and Z coordinates of the origin

          r
            Radius of the circle

          cr, cg, cb
            Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].

          w
            Line width of the circle

    RETURNS
          the CGO object (it also loads it into PyMOL, too).

```

### circleSelection
```
    circleSelection -- draws a cgo circle around a given selection or object

    PARAMS
          selName
            Name of the thing to encircle.

          r
            Radius of circle.
            DEFAULT: This cript automatically defines the radius for you.  If
            you select one atom and the resultant circle is too small, then
            you can override the script's calculation of r and specify your own.

          cr, cg, cb
            red, green and blue coloring, each a value in the range [0.0, 1.0]

    RETURNS
          The circle object.

```
 

## cgo_arrow.py

### cgo_arrow
```
    DESCRIPTION

        Create a CGO arrow between two picked atoms.

    ARGUMENTS

        atom1 = string: single atom selection or list of 3 floats {default: pk1}

        atom2 = string: single atom selection or list of 3 floats {default: pk2}

        radius = float: arrow radius {default: 0.5}

        gap = float: gap between arrow tips and the two atoms {default: 0.0}

        hlength = float: length of head

        hradius = float: radius of head

        color = string: one or two color names {default: blue red}

        name = string: name of CGO object
```
 

## cgo_grid.py

### eval_color
```None```

### cgo_grid
```
    DESCRIPTION

        Generates an animated flowing mesh object using the points provided
        or the current view. The shape is affected substantially by the arguments!

    USEAGE

        cgo_grid [ pos1 [, pos2 [, pos3 [, length_x [, length_z
                 [, npoints_x [, npoints_z [, nwaves_x [, nwaves_z
                 [, offset_x [, offset_z [, gain_x [, gain_z [, thickness
                 [, color [, nstates [, startframe [, endframe [, mode
                 [, view [, name [, quiet ]]]]]]]]]]]]]]]]]]]]]]

    EXAMPLE

        cgo_grid view=1

    ARGUMENTS

        pos1 = single atom selection (='pk1') or list of 3 floats {default: [0,0,0]}

        pos2 = single atom selection (='pk2') or list of 3 floats {default: [1,0,0]}

        pos3 = single atom selection (='pk3') or list of 3 floats {default: [0,0,1]}

        --> the plane is defined by pos1 (origin) and vectors to pos2 and pos3, respectively

        length_x = <float>: length of membrane {default: 30}
        length_z = <float>: length of membrane {default: ''} # same as length_x

        npoints_x = <int>: number of points(lines) along x-direction
                    {default: ''} #will be set to give a ~1 unit grid
        npoints_z = <int>: number of points(lines) along z-direction
                    {default: ''} #will be set to give a ~1 unit grid
                    {minimum: 1 # automatic}

        nwaves_x =   <float>: number of complete sin waves along object x-axis
                     {default: 2}
        nwaves_z =  <float>: number of complete sin waves along object z-axis
                    {default: ''} # same as nwaves_x
                    define separately to adjust number of waves in each direction



        offset_x = <float> phase delay (in degrees) of sin wave in x-axis
                 can be set to affect shape and starting amplitude {default: 0}
        offset_z = <float> phase delay (in degrees) of sin wave in z-axis
                 can be set to affect shape and starting amplitude
                 {default: ''} # same as  offset_x
        offset_x and offset_z can be used together to phase
        otherwise identical objects

        gain_x = <float>: multiplication factor for y-amplitude for x-direction
                 {default: 1}
        gain_z = <float>: multiplication factor for y-amplitude for z-direction
                 {default: ''} #=gain_x

        thickness = <float>: line thickness {default: 2}

        color = color name <string> (e.g. 'skyblue') OR
                rgb-value list of 3 floats (e.g. [1.0,1.0,1.0]) OR
                {default: ''} // opposite of background
                input illegal values for random coloring

        nstates =  <int>: number of states; {default: 60}
                   this setting will define how many states
                   the object will have (per wave) and how fluent and fast the
                   animation will be.
                   Higher values will promote 'fluent' transitions,
                   but decrease flow speed.
                       Note: Frame animation cycles thought the states one at a time
                       and needs to be set accordingly. Can also be used to phase
                       otherwise identical objects.
                   Set to 1 for static object {automatic minimum}

        startframe: specify starting frame <int> or set (='') to use current frame
                    set to 'append' to extend movie from the last frame {default: 1}
          endframe: specify end frame <int> or set (='') to use last frame
                    if 'append' is used for startframe,
                    endframe becomes the number of frames to be appended instead
                    {default: 1}
                    Note: if start- and endframe are the same, movie animation will
                    be skipped, the object will be loaded and can be used afterwards

        mode: defines positioning {default: 0}:
        0: pos1 is center
        1: pos1 is corner

        view {default: 0}:
        '0': off/ uses provided points to create CGO
        '1': overrides atom selections and uses current orienatation for positioning
             - pos1 = origin/center
             - pos2 = origin +1 in camera y
             - pos3 = origin +1 in camera z

        name: <string> name of cgo object {default: ''} / automatic

        quiet: <boolean> toggles output

```
 

## color_by_conservation.py
 

## colorblindfriendly.py
 

## colorbydisplacement.py

### ColorByDisplacementAll
```None```
 

## colorbyrmsd.py

### colorbyrmsd
```
    DESCRIPTION

        Align two structures and show the structural deviations in color to more
        easily see variable regions.

        Colors each mobile/target atom-pair by distance (the name is a bit
        misleading).

        Modifies the B-factor columns in your original structures.

    ARGUMENTS

        mobile = string: atom selection for mobile atoms

        target = string: atom selection for target atoms

        doAlign = 0 or 1: Superpose selections before calculating distances
        {default: 1}

        doPretty = 0 or 1: Show nice representation and colors {default: 1}

    EXAMPLE

        fetch 1ake 4ake, async=0
        remove chain B
        colorbyrmsd 1ake, 4ake
```
 

## cubes.py

### cubes
```
    DESCRIPTION

        Create a cube representation CGO for all atoms in selection.

    ARGUMENTS

        selection = string: atom selection {default: all}

        name = string: name of CGO object to create

        state = int: object state {default: 0 = all states}

        scale = float: scaling factor. If scale=1.0, the corners of the cube will
        be on the VDW surface of the atom {default: 0.5}

        atomcolors = 0/1: use atom colors (cannot be changed), otherwise
        apply one color to the object (can be changed with color command)
        {default: 1}

    SEE ALSO

        tetrahedra
```

### tetrahedra
```
    DESCRIPTION

        Create a tetrahedra representation CGO for all atoms in selection.

    SEE ALSO

        cubes
```
 

## cyspka.py

### cyspka
```None```

### loopcyspka
```None```
 

## displacementmap.py

### dispmap
```None```

### Coord
```None```
 

## distancetoatom.py

### distancetoatom
```
    DESCRIPTION

        distancetoatom.py
        Described at: http://www.pymolwiki.org/Distancetoatom

        Prints all distanced between the specified atom/coordinate/center
        and all atoms within cutoff distance that are part of the selection.
        All coordinates and distances can be saved in a csv-style text file report
        and can be appended to a (custom) atom property, if defined.

    USAGE

        distancetoatom [ origin [, cutoff [, filename [, selection
        [, state [, property_name [, coordinates [, decimals [, sort
        [, quiet ]]]]]]]]]]

    ARGUMENTS

        NAME        TYPE    FUNCTION
        origin:     <list>  defines the coordinates for the origin and can be:
                    <str>   1. a list with coordinates [x,y,z]
                            2. a single atom selection string {default='pk1'}
                            3. a multi-atom selection string (center will be used)
        cutoff      <float> sets the maximum distance {default: 10}
        filename    <str>   filename for optional output report. {default=None}
                            set to e.g. 'report.txt' to create a report
                            (omit or set to '', None, 0 or False to disable)
        selection   <str>   can be used to define/limit the measurment to specific
                            sub-selections {default='all'}
        state       <int>   object state, {default=0} # = current
        property_name <str> the distance will be stored in this property {p.dist}
                            set "" to disable
        coordinates <int>   toggle whether atom coordinated will be reported {0}
        decimals    <int>   decimals for coordinates and distance:
                            default = 3 # = max. PDB resolution
        sort        <int>   Sorting by distance?
                             1: ascending (default)
                             0: no sorting (by names)
                            -1: descending
        quiet       <bool>  toggle verbosity
```
 

## draw_rotation_axis.py
 

## drawgridbox.py

### drawgridbox
```
    DESCRIPTION
        Given selection, draw a grid box around it.

    USAGE:
        drawgridbox [selection, [nx, [ny, [nz, [padding, [lw, [r, [g, b]]]]]]]]

    PARAMETERS:
        selection,    the selection to enboxen
                      defaults to (all)

        nx,           number of grids on axis X
                      defaults to 10

        ny,           number of grids on axis Y
                      defaults to 10

        nz,           number of grids on axis Z
                      defaults to 10

        padding,      defaults to 0

        lw,           line width
                      defaults to 2.0

        r,            red color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        g,            green color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        b,            blue color component, valid range is [0.0, 1.0]
                      defaults to 1.0

    RETURNS
        string, the name of the CGO box

    NOTES
        * This function creates a randomly named CGO grid box. The user can
        specify the number of grids on X/Y/Z axis, the width of the lines,
        the padding and also the color.
```
 

## dssr_block.py

### dssr_block
```
    DESCRIPTION

        Create a nucleic acid base "block" cartoon with DSSR.

        Requires the "x3dna-dssr" program, available from http://x3dna.org/

    USAGE

        dssr_block [ selection [, state [, block_file [, block_depth
            [, block_color [, name [, exe ]]]]]]]

    ARGUMENTS

        selection = str: atom selection {default: all}

        state = int: object state (0 for all states) {default: -1, current state}

        block_file = face|edge|wc|equal|minor|gray: Corresponds to the --block-file
        option (see DSSR manual). Values can be combined, e.g. "wc-minor".
        {default: face}

        block_depth = float: thickness of rectangular blocks {default: 0.5}

        block_color = str: Corresponds to the --block-color option (new in DSSR
        v1.5.2) {default: }

        name = str: name of new CGO object {default: dssr_block##}

        exe = str: path to "x3dna-dssr" executable {default: x3dna-dssr}

    EXAMPLE

        fetch 1ehz, async=0
        as cartoon
        dssr_block
        set cartoon_ladder_radius, 0.1
        set cartoon_ladder_color, gray
        set cartoon_nucleic_acid_mode, 1

        # multi-state
        fetch 2n2d, async=0
        dssr_block 2n2d, 0
        set all_states

        # custom coloring
        fetch 1msy, async=0
        dssr_block block_color=N red | minor 0.9 | major yellow
```
 

## dynamic_mesh.py

### dynamic_mesh
```
    DESCRIPTION

        Make 'dynamic' mesh from volumetric data such as electron density map.
        The mesh will dynamically follow the center of the view.
        Contour level of isomesh can be changed by PageDown and PageUp keys.

        NOTE: Crystallographic operations can be applied to the map.

    USAGE

        dynamic_mesh map_name [, level [, radius [, name [, sym_source ]]]]

    ARGUMENTS

        map_name = string: name of volumetric object(map) to display

        level = float: contour level of isomesh {default: 1.0}

        radius = float: radius of isomesh around the center of the view {default: 8}

        name = string: name of mesh object {default: dynamic_mesh}

        sym_source = string: name of object from which symmetry
                             information is derived {default: map_name}

    EXAMPLE

        fetch 1hwk, async=0
        fetch 1hwk, 1hwk_map, type=2fofc, async=0
        dynamic_mesh 1hwk_map

    SEE ALSO

        isomesh
```
 

## elbow_angle.py

### elbow_angle
```

    DESCRIPTION

        Calculates the integer elbow angle of an antibody Fab complex and
        optionally draws a graphical representation of the vectors used to
        determine the angle.

    ARGUMENTS

        obj = string: object

        light/heavy = strings: chain ID of light and heavy chains, respectively

        limit_l/limit_h = integers: residue numbers of the last residue in the
        light and heavy chain variable domains, respectively

        draw = boolean: Choose whether or not to draw the angle visualization

    REQUIRES: com.py, transformations.py, numpy (see above)


```
 

## ex.py

### ex
```None```
 

## extra_fit.py

### extra_fit
```
    DESCRIPTION

        Like "intra_fit", but for multiple objects instead of
        multiple states.

    ARGUMENTS

        selection = string: atom selection of multiple objects {default: all}

        reference = string: reference object name {default: first object in selection}

        method = string: alignment method (command that takes "mobile" and "target"
        arguments, like "align", "super", "cealign" {default: align}

        ... extra arguments are passed to "method"

    SEE ALSO

        alignto, cmd.util.mass_align, align_all.py from Robert Campbell
```
 

## findSurfaceCharge.py

### findSurfaceAtoms
```
        Adapted from Jason Vertrees https://pymolwiki.org/index.php/FindSurfaceResidues
    DESCRIPTION

        Finds those atoms on the surface of a protein
        that have at least 'cutoff' exposed A**2 surface area.

    USAGE

        findSurfaceAtoms [ selection, [ cutoff ]]

    SEE ALSO

        findSurfaceResidues
```

### findSurfaceCharge
```
    DESCRIPTION

        Calculates a surface charge at entered pH. Also allows for the charge of an unfolded protein to be calculated.

    USAGE

        findSurfaceCharge [pH, [folded, [selection ,[cutoff]]]]

    ARGUMENTS

        pH = The pH value to estimate a surface charge at

        folded = Whether the protein is folded (True) or denatured (False)

        selection = string: object or selection in which to find exposed
            residues {default: empty string - all objects}

        cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

    RETURNS

        A printout of the estimated surface charge at a given pH

```
 

## findSurfaceResidues.py

### findSurfaceAtoms
```
    DESCRIPTION

        Finds those atoms on the surface of a protein
        that have at least 'cutoff' exposed A**2 surface area.

    USAGE

        findSurfaceAtoms [ selection, [ cutoff ]]

    SEE ALSO

        findSurfaceResidues
```

### findSurfaceResidues
```
    DESCRIPTION

        Finds those residues on the surface of a protein
        that have at least 'cutoff' exposed A**2 surface area.

    USAGE

        findSurfaceResidues [ selection, [ cutoff, [ doShow ]]]

    ARGUMENTS

        selection = string: object or selection in which to find exposed
        residues {default: all}

        cutoff = float: cutoff of what is exposed or not {default: 2.5 Ang**2}

    RETURNS

        (list: (chain, resv ) )
            A Python list of residue numbers corresponding
            to those residues w/more exposure than the cutoff.

```
 

## findseq.py

### findseq
```
    DESCRIPTION:
        Given a sequence/regex to find, select those matching amino acids in the
        protein.

    USAGE:
        findseq needle, [haystack[, selName[, het[, matchMode]]]]

    PARAMS:
    needle (string)
            the sequence of amino acids to match and select
            in the haystack.  This can be a sequence of amino
            acids, or a string-style regular expression.  See
            examples.

    haystack (string; PyMOL object or selection; defaults to *)
            name of the PyMOL object/selection in which
            to find the needle.

    selName (string; defaults to None)
            This is the name of the selection to return.  If selName
            is left blank (None), then the selection name will be
            foundSeqXYZ where XYZ is some random number; if selName is
            "sele" the usual PyMOL "(sele)" will be used; and, lastly,
            if selName is anything else, that name will be used verbatim.

    het (0 or 1; defaults to 0)
            This boolean flag allows (1) or disallows (0) heteroatoms
            from being considered.

    matchMode (first/all/chain; defaults to "all")
            Subsequences or motifs might be repeated, this controls how we
            consider multiple matches. Options are:
            - 'first': Return the first match found in each object.
            - 'chain': Return the first match found in each chain.
            - 'all': Return all matches found in all chains.

    RETURNS:
        a newly created selection with the atoms you sought.

    EXAMPLE:
        # find SPVI in 1h12, foundSeqXYZ as return name
        findseq SPVI, 1h12

        # find FATEW and make it (sele).
        findseq FATEW, 1g01, sele

        # find the regular expression GMS.*QWY in 1a3h
        # and put the return value in (sele).
        fetch 1a3h
        # this ends up finding the sequence, GMSSHGLQWY
        findseq GMS.*QWY, 1a3h, sele
```
 

## flatten_obj.py

### flatten_obj
```
    DESCRIPTION

        "flatten_obj" combines multiple objects or states into a single object,
        renaming chains where required

    USAGE

        flatten_obj name, selection[, state[, rename[, quiet[, chain_map]]]]

    ARGUMENTS

        name = a unique name for the flattened object {default: flat}

        selection = the set of objects to include in the flattening. The selection
            will be expanded to include all atoms of objects. {default: all}

        state = the source state to select. Use 0 or -1 to flatten all states {default: 0}

        rename = The scheme to use for renaming chains: {default: 0}
            (0) preserve chains IDs where possible, rename other chains
                alphabetically
            (1) rename all chains alphabetically
            (2) rename chains using the original chain letter, object name, and state

        quiet = If set to 0, print some additional information about progress and
            chain renaming {default: 1}

        chain_map = An attribute name for the 'stored' scratch object. If
            specified, `stored.<chain_map>` will be populated with a dictionary
            mapping the new chain names to a tuple giving the originated object,
            state, and chainID. {default: ""}

    NOTES

        Like the select command, if name is omitted then the default object name
        ("flat") is used as the name argument.

        Chain renaming is tricky. PDB files originally limited chains to single
        letter identifiers containing [A-Za-z0-9]. When this was found to be
        limiting, multi-letter chains (ideally < 4 chars) were allowed. This is
        supported as of PyMOL 1.7. Earlier versions do not accept rename=2, and
        will raise an exception when flattening a structure with more than 62
        chains.

    EXAMPLES

        flatten_obj flat, nmrObj
        flatten_obj ( obj1 or obj2 )

    SEE ALSO

        split_states

```
 

## focal_blur.py

### FocalBlur
```
    DESCRIPTION

        Creates fancy figures by introducing a focal blur to the image. The object
        at the origin will be in focus.

    AUTHOR

        Jarl Underhaug
        University of Bergen
        jarl_dot_underhaug_at_gmail_dot_com

        Updates by Jason Vertrees and Thomas Holder

    USAGE

        FocalBlur aperture=float, samples=int, ray=0/1, width=int, height=int

    EXAMPELS

        FocalBlur aperture=1, samples=100
        FocalBlur aperture=2, samples=100, ray=1, width=600, height=400
```
 

## format_bonds.py

### format_bonds
```
    DESCRIPTION
        Formats bonds in aromatic or charged residues
    EXAMPLE
        frag PHE
        format_bonds
    USAGE
        format_bonds [ selection [, bonds ]]
    ARGUMENTS
        selection: <str> input selection {default: 'all'}
        bonds:     <int> toogles format of bonds
                   1: single bonds (deactivates valence display)
                   2: regular double bonds (activates valence display)
                 >=3: delocalized (activates valence display)
```
 

## forster_distance_calculator.py

### forster
```None```
 

## get_colors.py

### get_colors
```
    DESCRIPTION:
        returns a list of available pymol colors

    USAGE:
        get_colors [ selection [, quiet ]]

    EXAMPLES:
        get_colors # basic colors
        get colors all # larger range with intermediates
```

### get_random_color
```
    DESCRIPTION:
        returns a random color name available in pymol
        ! Requires get_colors !Indended mostly for use in Python

    USAGE:
        get_random_color [ selection [, quiet ]]

    EXAMPLES:
        # print a random color name:
        get_random_color
        # color object randomly:
        fetch 1hpv, async=0
        cmd.color(get_random_color())
```
 

## get_raw_distances.py

### get_raw_distances
```
    DESCRIPTION

        Get the list of pair items from distance objects. Each list item is a
        tuple of (index1, index2, distance).

        Based on a script from Takanori Nakane, posted on pymol-users mailing list.
        http://www.mail-archive.com/pymol-users@lists.sourceforge.net/msg10143.html

    ARGUMENTS

        names = string: names of distance objects (no wildcards!) {default: all
        measurement objects}

        state = integer: object state {default: 1}

        selection = string: atom selection {default: all}

    SEE ALSO

        select_distances, cmd.find_pairs, cmd.get_raw_alignment
```

### select_distances
```
    DESCRIPTION

        Turns a distance object into a named atom selection.

    ARGUMENTS

        names = string: names of distance objects (no wildcards!) {default: all
        measurement objects}

        name = a unique name for the selection {default: sele}

    SEE ALSO

        get_raw_distances
```
 

## grepset.py

### grepset
```
    DESCRIPTION
        "grepset" greps through the list of settings using a python
        regular expression as defined in the 're' module.
        It returns a list of settings/values matching the regexp.
        No regexp returns every setting.

    USAGE
        grepset [regexp]

    EXAMPLE
        grepset line
        grepset ray
        grepset (^line|color$)

    SEE ALSO
            Python re module
```
 

## gridbox.py

### gridbox
```
    DESCRIPTION
    Create a box from the center coordinate of the box and the size of box

    USAGE
    run gridbox.py
    1the simplest
    gridbox center_x,center_y,center_z,size_x,size_y,size_z
    2rename the box object
    gridbox center_x,center_y,center_z,size_x,size_y,size_z,name,
    3set the color of the box object
    gridbox center_x,center_y,center_z,size_x,size_y,size_z,name,r1,g1,b1
    4set the trasp of the box
    gridbox center_x,center_y,center_z,size_x,size_y,size_z,name,r1,g1,b1,trasp

    ps:the value of r1,g1,b1 trasp   range  is 0-1
       trasp=1,no trasprent


```
 

## hbplus.py

### hbplus
```
    DESCRIPTION

        HBPLUS wrapper
```
 

## inertia_tensor.py

### tensor
```
    DESCRIPTION

        This script will draw the inertia tensor of the selection.

    ARGUMENTS

        selection = string: selection for the atoms included in the tensor calculation

        name = string: name of the tensor object to be created {default: "tensor"}

        state = int: state/model in the molecule object used in the tensor calculation

        scaling = int {0, 1, or 2}: 0 for no scaling of the inertia axes, 1 for scaling
        according to the molecular shape, 2 for scaling according to the eigenvalues 
        {default: 0}

    EXAMPLE

        PyMOL> run inertia_tensor.py
        PyMOL> tensor molecule_object & i. 2-58+63-120 & n. n+ca+c, "tensor_model5_dom2", 5, 1

    NOTES

        Requires numpy.
```
 

## isoslider.py

### isoslider
```
    DESCRIPTION

        Opens a dialog with isolevel sliders for all isomesh and isosurface
        objects in PyMOL.
```
 

## load_img_stack.py

### load_img_stack
```
    DESCRIPTION

        Load a stack of images as a map

    ARGUMENTS

        pattern = str: image filename or pattern

        name = str: map object name to create

        grid = float: grid spacing in Angstrom {default: 1.0}

        channel = int: color channel for RGB images {default: 0}

        normalize = 0 or 1: normalize data {default: 1}

        extent = 3-float: (a,b,c) edge lengths in Angstrom, overwrites "grid"
        arguments if given {default: }

    EXAMPLES

        load_img_stack img*.tif, extent=(21.0, 14.5, 18.2)
```
 

## load_nrrd.py
 

## modevectors.py

### modevectors
```
    Authors Sean Law & Srinivasa
    Michigan State University
    slaw_(at)_msu_dot_edu

    Editor Sacha Yee

    USAGE

    While in PyMOL

    Parameter                Preset            Type    Description
    first_obj_frame          Undefined         String  Object name of the first structure.  The mode vector will propagate from this structure.  Defined by user.
    last_obj_frame           Undefined         String  Object name of the last structure.  The mode vector (arrow head) will end at this structure.  Defined by user.
    first_state              1                 Integer Defines state of first object
    last_state               1                 Integer Defines state of last object
    outname                  modevectors       String  Name of object to store mode vectors in.
    head                     1.0               Float   Radius for the circular base of the arrow head (cone)
    tail                     0.3               Float   Radius for the cylinder of the arrow tail (cylinder)
    head_length              1.5               Float   Length of the arrow head (from the base of the cone to the tip of cone)
    head_rgb                 1.0,1.0,1.0       String  RGB colour for the arrow head.
    tail_rgb                 1.0,1.0,1.0       String  RGB colour for the arrow tail.
    cutoff                   4.0               Float   Skips mode vectors that do not meet the cutoff distance (in Angstroms).
    skip                     0                 Integer Denotes how many atoms to skip.  No arrows will be created for skipped atoms.
    cut                      0.0               Float   Truncates all arrow tail lengths (without disturbing the arrow head) (in Angstroms).
    atom                     CA                String  Designates the atom to derive mode vectors from.
    stat                     show              String  Keeps track and prints statistics (total modevectors, skipped, cutoff).
    factor                   1.0               Float   Multiplies each mode vector length by a specified factor.
                                                       Values between 0 and 1 will decrease the relative mode vector length.
                                                       Values greater than 1 will increase the relative mode vector length.
    notail                   0                 Integer Hides tails and only uses cones (porcupine plot)
```
 

## movie_color_fade.py

### movie_color_fade
```
    DESCRIPTION

        Fades the color of representations in movies
        #NB!: Defines and uses new color names using the selection name and frame numbers

    USE

        movie_color_fade startframe='', startcolor=red, endframe='', endcolor=green, selection=all
        #defaults indicated

    PARAMETERS

       startframe, endframe = beginning and end movie frame for fading
       startcolor, endcolor = coloring at start and end
       selection: target selection

       NB! startframe and endframe can be omitted or set='' to assign current and last frame respectively

    EXAMPLES

        ##### 1. #####
        fetch 1hpv, async=0
        as cartoon
        orient
        color yellow
        mset 1x120
        movie_color_fade 1, yellow, 60, blue
        movie_color_fade 60, blue, 120, yellow
        #####

        ##### 2. #####
        #repeat command and specify 'selection' to change multiple colors
        fetch 1hpv, async=0
        as cartoon
        orient
        color white
        mset 1x60
        movie_color_fade auto,white,auto,skyblue,ss s
        movie_color_fade auto,white,auto,red,ss h
        movie_color_fade auto,white,auto,grey,ss l+""
        #####

    SEE ALSO

        mdo, mappend, set, movie_fade
```
 

## movie_fade.py

### movie_fade
```
    DESCRIPTION

        Fades representations in movies with their transparency settings.

    USAGE

        movie_fade setting, startFrame, startVal, endFrame, endVal [, selection ]

    EXAMPLE

        fetch 1rx1, async=0
        as cartoon
        show surface
        mset 1x80
        movie.roll
        movie_fade transparency,  1, 0., 40, 1.
        movie_fade transparency, 41, 1., 80, 0.

    SEE ALSO

        mdo, mappend, set
```
 

## nmr_cnstr.py
 

## perp_maker.py

### perp_maker
```
    DESCRIPTION

        Creates perpendicular planes
```
 

## plane.py

### make_plane
```
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, a1, a2, a3

    where each atom is a standard PyMOL selection (defaults to pk1,pk2 and pk3)
```

### make_plane_points
```
    DESCRIPTION
    Create a CGO plane from three atomic coordinates

    USAGE
    make_plane name, l1, l2, l3

    where each xys is a list with floats of x,y,z coordinates
```
 

## plot_noe.py

### plot_noe
```
    DESCRIPTION

        A function for plotting XPLOR NOE restraints on a structure

    ARGUMENTS

        filename = string: The filename of the NOE retraint file in XPLOR NIH format.

        line_color = string: The color for the NOE lines. {default: yellow}

        line_width = float: The thickness of the NOE lines. {default: 1.0}

        advanced_coloring = color restraints by distance.

        single = string: create a single object for all restraints.

        aria = integer: Name NOEs after Aria IDs.

        per_atom: Group NOEs on atom basis

        per_residue: Group NOEs on residue basis (default)

    NOE Restraint Format Example

        assign (residue 5 and name HB#) (residue 21 and name HA) 3.0 0.7 0.7

    EXAMPLE

        PyMOL> plot_noe noe_short.tbl
```
 

## poseview.py

### poseview
```
    DESCRIPTION

        PoseView wrapper

        http://www.biosolveit.de/poseview/

    USAGE

        poseview [ ligand [, protein [, width [, height [, exe [, state ]]]]]]

    ARGUMENTS

        ligand = string: atom selection {default: organic inorganic}

        protein = string: atom selection {default: polymer}

        width = int: image width {default: viewport width}

        height = int: image height {default: viewport height}

        filename = string: output PNG file name {default: temporary}

        exe = string: path to executable {default: poseview}

    SETUP

        1) Put poseview executable to PATH (e.g. /usr/bin/poseview)
        2) Set environment variable BIOSOLVE_LICENSE_FILE="/path/to/poseview.lic"
```
 

## propka.py

### propka
```    20110823```

### getpropka
```None```
 

## pymol2glmol.py

### pymol2glmol
```None```
 

## pymolscriptrepo.py
 

## quickdisplays.py

### disp_list
```None```

### disp_ss
```
    DESCRIPTION

        Formats the passed object into secondary structure cartoon

    USAGE

        disp_ss [ selection [, colors [, only ]]]

    PARAMETERS

        NAME=DEFAULT               TYPE    FUNCTION
        selection='all'            <str>   input selection
        colors='marine red white'  <str>   any three colors for: sheets, helices and loops
                                           e.g. 'marine red white'
                                           can also be set to either util.cbc, util.rainbow,
                                           or util.chainbow (alone)
                                           set='False' to supress coloring altogether, or enter False
                                           for the coloring to be omitted, e.g. 'marine False green'
        only                       <bool>  if True will use show_as; else show
```

### disp_ball_stick
```
    DESCRIPTION

        Formats the passed object into ball and stick

    USEAGE

        disp_ball_stick [ selection [, hydrogens [, only ]]]

    EXAMPLE

        fetch 1hpv, async=0
        disp_ball_stick
        util.cbaw

    PARAMETERS

        NAME=DEFAULT       TYPE    FUNCTION
        selection='all'    <str>   input selection
        hydrogens          <int>   -1: remove; 1: add; else: as is
        only=False         <bool>  if True will use show_as; else show

```

### disp_stick_ball
```
    see help disp_stick_ball
```

### disp_mesh
```
    DESCRIPTION

        Adds a mesh to the object
        Has advanced coloring options and automatically accounts for the hydrogens

    USEAGE

        disp_mesh [ selection [, color_m [, hydrogens [, only [, limits]]]]]
        disp_mesh selection=all, color_m=default
        disp_mesh selection=all, color_m=white
        disp_mesh selection=all, color_m=putty

    PARAMETERS

        NAME=DEFAULT       TYPE    FUNCTION
        selection='all'    <str>   input selection
        color_m='default'  <str>   'default': as current
                                   'name': colors by color or ramp called name
                                   'putty': b-factor on surface
        hydrogens=0        <int>   -1: remove; 1: add; else: as is
        only=False         <bool>  if True will use show_as; else show
        limits=5           <list or flaot>
                                   applies only if color_m=='putty'
                                   sets the b-factor range limits
                                   <list> [min,max] # absolute values
                                   <float> percentile cutoff (both sides) # relative for each protein
```

### disp_surf
```
    DESCRIPTION

        Advanced surface representation (cf. examples)

    USAGE

        disp_surf [ selection [, color_s [, transparency [, hydrogens [, solvent [, ramp_above [, only [, limits]]]]]]]]

    EXAMPLES

        disp_surf # opaque surface with default colors
        disp_surf all, white, 0.5 # half-transparent white surface
        disp_surf all, putty # b-factor on surface

    PARAMETERS

        NAME=DEFAULT       TYPE    FUNCTION
        selection='all'    <str>   input selection
        color_s='default'  <str>   'default': as current
                                   'name': colors by color or ramp called name
                                   'putty': b-factor on surface (by resi)
        transparency=0     <float> set surface transparency
        hydrogens=0        <int>   -1: remove; 1: add; else: as is
        solvent=0          <int>   defines 'surface_solvent'
        ramp_above=1       <int>   defines 'surface_ramp_above_mode'
        only=False         <bool>  if True will use show_as; else show
        limits=5           <list or flaot>
                                   applies only if color_s=='putty'
                                   sets the b-factor range limits
                                   <list> [min,max] # absolute values
                                   <float> percentile cutoff (both sides) # relative for each protein
```

### disp_putty
```
    DESCRIPTION

        Formats the passed object into a Putty b-factor sausage

    USEAGE

        disp_putty [ selection ]
        selection    <str>    input selection
        limits=10    <list or flaot>
                              applies only if color_m=='putty'
                              sets the b-factor range limits (by protein)
                              <list> [min,max]
                              <float> percentile cutoff (both sides)
        only=True             <bool>  if True will use show_as; else show
```
 

## removealt.py

### removealt
```
    removeAlt -- remove all alternate location-atoms not of altloc "keep" from object.

    input:
            obj -- the object(s) to remove the atoms frmo
            keep -- which type of alt loc to keep

    output: none -- removes atoms

    examples:
            removeAlt # remove all altLocations that aren't altloc A
            removeAlt pdbID, C  # remove all but C altlocations from pdbID
```
 

## renumber.py

### renumber
```
    DESCRIPTION

        Set residue numbering (resi) based on connectivity.

    ARGUMENTS

        selection = string: atom selection to renumber {default: all}

        start = integer: counting start {default: 1}

        startsele = string: residue to start counting from {default: first in
        selection}
```
 

## rotkit.py

### rotateline
```None```

### mutate
```None```

### toline
```None```
 

## save_mopac.py

### save_mopac
```
    DESCRIPTION

        Save to MOPAC format

    ARGUMENTS

        filename = string: file path to be written

        selection = string: atoms to save {default: all}

        zero = string: atoms to save with zero flag {default: none}

        state = integer: state to save {default: -1 (current state)}
```
 

## save_pdb_with_anisou.py

### save_pdb_with_anisou
```
    DESCRIPTION

         Save in PDB format including ANISOU records.

    SEE ALSO

        save
```
 

## save_settings.py

### save_settings
```
    DESCRIPTION

        Dumps all settings with non-default values to ~/.pymolrc-settings.py

        Feature Request: Save settings for later use - ID: 1009951
        https://sourceforge.net/tracker/?func=detail&aid=1009951&group_id=4546&atid=354546
```
 

## select_sites.py

### select_sites
```
    DESCRIPTION

        Make named selections from SITE records.

    ARGUMENTS

        name = string: molecular object {default: all}

        filename = string: PDB file name with SITE records {default: look in
        current directory and fetch_path for <name>.pdb}

        prefix = string: prefix for named selections {default: site_}

        nice = 0 or 1: make colored sticks representation for sites {default :1}
```

### sites
```None```
 

## show_bumps.py

### show_bumps
```
    DESCRIPTION

        Visualize VDW clashes

    ARGUMENTS

        selection = string: atom selection {default: all}

        name = string: name of CGO object to create {default: bump_check}
```
 

## show_ligand_interactions.py

### show_ligand_interactions
```
    DESCRIPTION

        Visualize interactions between receptor and ligand.

    ARGUMENTS

        recsel = string: atom selection of the receptor {default: "not hetatm"}

        ligsel = string: atom selections of the ligand {default: "hetatm"}

        cutoff = float: show as sticks all receptor residues within this distance from the ligand {default: 5.0}
```
 

## spectrum_states.py

### spectrum_states
```
    DESCRIPTION

        Color each state in a multi-state object different.

    USAGE

        spectrum_states [ selection [, representations [, color_list [, first [, last ]]]]]

    ARGUMENTS

        selection = string: object names (works with complete objects only)
        {default: all}

        representations = string: space separated list of representations
        {default: cartoon ribbon}

        color_list = string: space separated list of colors {default: blue cyan
        green yellow orange red}

    SEE ALSO

        spectrum, spectrumany
```
 

## spectrumany.py

### spectrumany
```
    DESCRIPTION

        Define a color spectrum with as many color-stops as you like (at least 2).

    USAGE

        spectrumany expression, color_list [, selection [, minimum [, maximum ]]]

    ARGUMENTS

        expression = count, resi, b, q, or pc: respectively, atom count, residue
        index, temperature factor, occupancy, or partial charge {default: count}

        color_list = string: Space separated list of colors

        ... all other arguments like with `spectrum` command

    EXAMPLE

        spectrumany count, forest green yellow white
        spectrumany b, red yellow white, (polymer), maximum=100.0

    SEE ALSO

        spectrum
```
 

## spectrumbar.py

### spectrumbar
```
    Author Sean M. Law
    University of Michigan
    seanlaw_(at)_umich_dot_edu

    USAGE

    While in PyMOL

    run spectrumbar.py

    spectrumbar (RGB_Colors,radius=1.0,name=spectrumbar,head=(0.0,0.0,0.0),tail=(10.0,0.0,0.0),length=10.0, ends=square)

    Parameter     Preset         Type     Description
    RGB_Colors    [1.0,1.0,1.0]  N/A      RGB colors can be specified as a
                                          triplet RGB value or as PyMOL
                                          internal color name (i.e. red)
    radius        1.0            float    Radius of cylindrical spectrum bar
    name          spectrumbar    string   CGO object name for spectrum bar
    head          (0.0,0.0,0.0)  float    Starting coordinate for spectrum bar
    tail          (10.0,0.0,0.0) float    Ending coordinate for spectrum bar
    length        10.0           float    Length of spectrum bar
    ends          square         string   For rounded ends use ends=rounded

    Examples:

    spectrumbar red, green, blue
    spectrumbar 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0

    The above two examples produce the same spectrumbar!

    spectrumbar radius=5.0
    spectrumbar length=20.0

```
 

## stereo_ray.py

### stereo_ray
```
    DESCRIPTION

       "stereo_ray" ray-traces the current scene twice (separated by 
       a six-degree rotation around the y axis)
       and saves a pair of images that can be combined in any image
       manipulation software to form a stereoimage.
       The first argument, the output file name, is mandatory.
       The second and third arguments, the size of the image, are not.
       If the width is given, the height will be calculated.

    USAGE

       stereo_ray filename [, width [, height]]

    EXAMPLE

       stereo_ray output, 1000, 600
       stereo_ray secondImage.png
```
 

## tmalign.py

### alignwithanymethod
```
    DESCRIPTION

        Align copies of mobile to target with several alignment methods

    ARGUMENTS

        mobile = string: atom selection

        target = string: atom selection

        methods = string: space separated list of PyMOL commands which take
        arguments "mobile" and "target" (in any order) {default: align super
        cealign tmalign}
```

### tmalign
```
    DESCRIPTION

        TMalign wrapper

        Reference: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9
        http://zhanglab.ccmb.med.umich.edu/TM-align/

    USAGE

        tmalign mobile, target [, args [, exe ]]

    ARGUMENTS

        mobile, target = string: atom selections

        args = string: Extra arguments like -d0 5 -L 100

        exe = string: Path to TMalign executable {default: TMalign}

        ter = 0/1: If ter=0, then ignore chain breaks because TMalign will stop
        at first TER record {default: 0}

    SEE ALSO

        tmscore, mmalign
```

### tmscore
```
    DESCRIPTION

        TMscore wrapper

        Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710
        http://zhanglab.ccmb.med.umich.edu/TM-score/

    ARGUMENTS

        mobile, target = string: atom selections

        args = string: Extra arguments like -d 5

        exe = string: Path to TMscore executable {default: TMscore}

        ter = 0/1: If ter=0, then ignore chain breaks because TMscore will stop
        at first TER record {default: 0}

    SEE ALSO

        tmalign, mmalign
```

### mmalign
```
    DESCRIPTION

        MMalign wrapper

        Reference: S. Mukherjee and Y. Zhang, Nucleic Acids Research 2009; 37: e83
        http://zhanglab.ccmb.med.umich.edu/MM-align/

    SEE ALSO

        tmalign, tmscore
```
 

## togroup.py

### toGroup
```
    DESCRIPTION
    toGroup will take a multistate object and extract it
    to a group with N objects all in state #1.  It essentially
    performs the following:

    split_states myObj, prefix=somePrefix
    group newGroup, somePrefix*
    delete myObj

    PARAMETERS:

    groupName (string)
        The name of the group to create

    sel (string)
        The name of the selection/object from which
        to make the group

    prefix (string)
        The prefix of the names of each of split states.
        For example, if your prefix is ''obj'' and is in
        states 1 through 100 then the states will be labeled
        obj1, obj2, obj3, ..., obj100.

    delOrig (string/boolean)
        If true then delete the original selection, otherwise not.

    RETURN

    Nothing, it makes a new group.

```
 

## transformations.py
 

## uniprot_features.py

### uniprot_features
```
    DESCRIPTION

        Fetch feature list from uniprot.org and create named selections.

        Requires residue numbering (resi) to match uniprot sequence!

    ARGUMENTS

        uniprot_id = string: UniProtKB name or accession

        selection = string: atom selection {default: all}

        withss = 0/1: update secondary structure {default: 0}
```

### uniprot_auto
```
    DESCRIPTION

        Like "uniprot_features" but with automatic fetching of UniProtKB accession
        and sequence mapping for given pdb_id from http://www.bioinf.org.uk/pdbsws/

    ARGUMENTS

        pdb_id = string: PDB accession ID

        selection = string: atom selection {default: <pdb_id>, will be fetched if
        no such object is loaded}

        withss = 0/1: update secondary structure {default: 0}
```
 

## viol_noes.py

### viol_noes
```
    DESCRIPTION

        Visualize Xplor-NIH NOE violations.

    ARGUMENTS

        molecule = string: molecule on which to show the violations.

        viol_file = string: Xplor-NIH .viol file that contains the violations to be visualized.

        viol_class = string: NOE class in .viol file to show {default: None (means all NOE classes)}.


    EXAMPLE

        PyMOL> run viol_noes.py
        PyMOL> viol_noes molecule, ./molecule.pdb.viol

    NOTES

        The NOE violations will be shown as distances between the relevant residues/atoms and colored according to the severity of violation (the closer to the blue end of the spectrum, the more severe the violation; to closer to the red end, the less severe the violation).
```
 

## wfmesh.py

### createWFObj
```None```
 

## gendocs.py
 

## plugins/dynoplot.py

### ramachandran
```None```
 

## plugins/get_box.py
 

## plugins/outline.py
 

## plugins/views.py
 

## plugins/vina.py
 
