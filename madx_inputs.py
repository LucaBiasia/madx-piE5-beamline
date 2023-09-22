end_sequence = '''
endsequence;
use, sequence=piE5;
'''

twiss = '''
select, flag=twiss;
use, sequence=piE5;
twiss, betx=sigmax^2/ex, bety=sigmay^2/ey, alfx=alfax, alfy=alfay, exact;

select, flag=aperture;
aperture;
// use sequence
use, sequence=piE5;
'''

ptc_twiss = f'''
use, sequence=piE5;

ptc_create_universe;
ptc_create_layout, time=false, model=2, method=4, nst=1, exact, closed_layout=false;
select, flag=ptc_twiss;
ptc_twiss, icase=4, no=5, savemaps=true, normal, betx=sigmax^2/ex, bety=sigmay^2/ey, alfx=alfax, alfy=alfay, closed_orbit=False;
ptc_select_moment, table=moments, moment_s=2,002,11,0011,02,0002,1,001;
ptc_moments, no=2, xdistr=gauss, ydistr=gauss;
ptc_end;
'''

ptc_create = '''
ptc_create_universe, ntpsa;
ptc_create_layout, time=false, model=2, method=4, nst=1, resplit, thin=0.001, closed_layout=false;    //time must be switched off!!!
ptc_setswitch, time=false, totalpath=false;
'''

ptc_track = '''
ptc_track, icase=5, element_by_element=true, turns=1, norm_no=5,
    dump=true, onetable=true, recloss=true, maxaper={0.425, 100, 0.425, 100, 100}, file=tracks_piE5_;

// ENDING TRACKING MODULE
ptc_track_end;
ptc_end;
'''