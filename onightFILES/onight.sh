#CALCULATE SOME MODELS OVER NIGHT
function run_model () {
    diff1d.eo < in > output_diff.log
    sc3d.eo < in > output_sc3d.log
}
#
function set_indat () {
    #calculate source functions
    output_file='indat_model.nml'
    eps_cont=$1
    xmin=$2
    xmax=$3
    ymin=$4
    ymax=$5
    opt_alo_cont=$6
    opt_method=$7
    echo $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
    
    echo '################# setting indat file####################'
    echo 'output_file: ' $output_file

    echo '&input_options' > $output_file
    echo 'opt_ng_cont=f' >> $output_file
    echo 'opt_ait_cont=f' >> $output_file
    echo 'opt_alo_cont='$opt_alo_cont >> $output_file
    echo 'opt_angint_method=0' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&input_model' >> $output_file
    echo 'rstar=10.d0' >> $output_file
    echo 'teff=40.d3' >> $output_file
    echo 'trad=40.d3' >> $output_file
    echo 'tmin=0.8d0' >> $output_file
    echo 'xmloss = 5.d-5' >> $output_file
    echo 'vmin = 10.d0' >> $output_file
    echo 'vmax = 2.d3' >> $output_file
    echo 'beta = 1.d0' >> $output_file
    echo 'yhe = 0.1d0' >> $output_file
    echo 'hei = 2.d0' >> $output_file
    echo 'xnue0=1.93798D15' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&input_cont' >> $output_file
    echo 'eps_cont ='$eps_cont >> $output_file
    echo 'kcont = 1.d0' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&dimensions_3dz' >> $output_file
    echo 'nz=31' >> $output_file
    echo 'zmin =' >> $output_file 0.d0
    echo 'zmax = 10.d0' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&dimensions_3dx' >> $output_file
    echo 'nx=31' >> $output_file
    echo 'xmin = '$xmin >> $output_file
    echo 'xmax = '$xmax >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&dimensions_3dy' >> $output_file
    echo 'ny=31' >> $output_file
    echo 'ymin = '$ymin >> $output_file
    echo 'ymax = '$ymax >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&dimensions_freq' >> $output_file
    echo 'nnue=1' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&dimensions_angles' >> $output_file
    echo 'ntheta=1' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&input_diff1d' >> $output_file
    echo 'opt_method=4' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&input_sc1d' >> $output_file
    echo 'opt_method=5' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&input_sc2d' >> $output_file
    echo 'opt_method=5' >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&input_sc3d' >> $output_file
    echo 'opt_method='$opt_method >> $output_file
    echo '/' >> $output_file
    echo '' >> $output_file
    echo '&benchmark' >> $output_file
    echo 'benchmark_mod = 0' >> $output_file
    echo 'theta = 85.' >> $output_file
    echo 'phi = 315. /' >> $output_file
    echo '' >> $output_file
    echo '' >> $output_file
    echo '&namelist_dum' >> $output_file
    echo '/' >> $output_file
    echo    
}
function copy_all () {

    modelstr=$1
    output_dir='/home/levin/Postdoc/notes/notes_bc_periodic/models3d/'$modelstr
    output_dir1=$output_dir'/sc3d'
    output_dir2=$output_dir'/diff1d'
    mkdir $output_dir
    mkdir $output_dir1
    mkdir $output_dir2
    cp outputFILES/sc3d/model3d.h5 $output_dir1
    cp output_sc3d.log $output_dir1
    cp outputFILES/diff1d/*.dat $output_dir2
    cp output_diff.log $output_dir2
}






#
#-------------------------------------------------------------------
#
modelstr='model0020'
eps_cont='1.d-5'
xmin='-5.d1'
xmax='5.d1'
ymin='-5.d1'
ymax='5.d1'
opt_alo_cont='3'
opt_method='4'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
#run_model
#copy_all $modelstr

exit 0







#
#-------------------------------------------------------------------
#
modelstr='model0300'
eps_cont='1.d-5'
xmin='-5.d1'
xmax='5.d1'
ymin='-5.d1'
ymax='5.d1'
opt_alo_cont='1'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0301'
eps_cont='1.d-5'
xmin='-5.d0'
xmax='5.d0'
ymin='-5.d0'
ymax='5.d0'
opt_alo_cont='1'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0302'
eps_cont='1.d-5'
xmin='-5.d-1'
xmax='5.d-1'
ymin='-5.d-1'
ymax='5.d-1'
opt_alo_cont='1'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0310'
eps_cont='1.d-5'
xmin='-5.d1'
xmax='5.d1'
ymin='-5.d1'
ymax='5.d1'
opt_alo_cont='2'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0311'
eps_cont='1.d-5'
xmin='-5.d0'
xmax='5.d0'
ymin='-5.d0'
ymax='5.d0'
opt_alo_cont='2'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0312'
eps_cont='1.d-5'
xmin='-5.d-1'
xmax='5.d-1'
ymin='-5.d-1'
ymax='5.d-1'
opt_alo_cont='2'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0320'
eps_cont='1.d-5'
xmin='-5.d1'
xmax='5.d1'
ymin='-5.d1'
ymax='5.d1'
opt_alo_cont='3'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0321'
eps_cont='1.d-5'
xmin='-5.d0'
xmax='5.d0'
ymin='-5.d0'
ymax='5.d0'
opt_alo_cont='3'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr
#
#-------------------------------------------------------------------
#
modelstr='model0322'
eps_cont='1.d-5'
xmin='-5.d-1'
xmax='5.d-1'
ymin='-5.d-1'
ymax='5.d-1'
opt_alo_cont='3'
opt_method='7'

set_indat $eps_cont $xmin $xmax $ymin $ymax $opt_alo_cont $opt_method
run_model
copy_all $modelstr



