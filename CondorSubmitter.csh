#! /bin/env csh

set Exec = ./bin/RunTest

# make sure executable exists
make $Exec || exit

####### Initialize condor file
echo ""  > CondorFile
echo "Universe     = vanilla" >> CondorFile
echo "Executable   = ${Exec}" >> CondorFile
echo "getenv = true" >>CondorFile
# echo "Notification = Complete" >> CondorFile
# echo "Notify_user  = kkauder@gmail.com"  >> CondorFile

set breit = 0

switch ($breit)
    case 0 : 
    set NameBase=Lab
    breaksw
    case 1 : 
    set NameBase=Breit
    breaksw
    default : 
    echo "unknown breit setting $breit"
    exit -1
endsw


# split into chunks
#set base = Data/extraSimu/PYTHIA/ep/TREES/pythia.ep.20x250.*Q2=10.0-100*root
set base = Data/extraSimu/PYTHIA/ep/TREES/pythia.ep.20x250.*root

foreach input ( ${base}* )
    # arguments
    set OutBase=`basename $input | sed 's/.root//g'`
    set OutName    = Results/Pieces/${NameBase}_${OutBase}.root
    set Files      = ${input}
    
    # Logfiles.
    set LogFile    = logs/${NameBase}_${OutBase}.out
    set ErrFile    = logs/${NameBase}_${OutBase}.err

    ### hand to condor
    set Args = ( -o $OutName -i $Files -breit $breit )
    echo "" >> CondorFile
    echo "Output       = ${LogFile}" >> CondorFile
    echo "Error        = ${ErrFile}" >> CondorFile
    echo "Arguments    = ${Args}" >> CondorFile
    echo "Queue" >> CondorFile   

    echo Submitting:
    echo $Exec $Args
    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile
    echo
end
condor_submit CondorFile


