#! /bin/bash
# 2017Feb16 install CDF patch into IDL and ENVI
VERSION="3.9.2"
ftp="https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf39_2/idl"
echo " "
echo "This is the script to install the latest IDL patch from SPDF onto your Linux/Unix"
echo "You may need to enter the administrator account's password."
echo "Only 64-bit of IDL is supported."
# <https://cdf.gsfc.nasa.gov/html/cdf_patch_for_idl.html>
IDL_Linux64="$ftp/linux/linux64/idl_cdf.so"
#IDL_Linux32="$ftp/linux/linux32/idl_cdf.so"
#IDL_Sparc64="$ftp/solaris/sparc64/idl_cdf.so"
#IDL_Sparc32="$ftp/solaris/sparc32/idl_cdf.so"
IDL_DLM="$ftp/idl_cdf.dlm"
LEAPS="https://cdf.gsfc.nasa.gov/html/CDFLeapSeconds.txt"
SHOWVERS="$ftp/showcdfversion.pro"
SHOWVERS2="$ftp/showcdfversion2.pro"
I32=""
if [ -x "$(command -v curl)" ]; then
  GET="curl -O "
elif [ -x "$(command -v wget)" ]; then
  GET="wget -c -nv --no-check-certificate "
else
  echo 'Error: require curl or wget' >&2
  exit 1;
fi
# find all installed versions of IDL and ENVI
echo " "
echo "This is the list of currently available IDL version(s) installed on"
echo "your Linux/Unix in /usr/local directory:"
echo " "
ToFind=""
if [ -d "/usr/local/nv5" ]; then
  ToFind="/usr/local/nv5"
fi
if [ -d "/usr/local/harris" ]; then
  ToFind="/usr/local/harris"
fi
if [ -d "/usr/local/exelis" ]; then
  ToFind="$ToFind /usr/local/exelis"
fi
if [ -d "/usr/local/itt" ]; then
  ToFind="$ToFind /usr/local/itt"
fi
options=( $(find ${ToFind} -name idl_cdf.so -exec dirname {} \;) )
tocont=0
while [[ ${tocont} == 0 ]] ; do
 PS3=$'\n'"Select the number, shown on the left of the above list, for the IDL/ENVI to be replaced by the CDF patch, or 0 to end: "
 select opt in "${options[@]}" ; do
  if (( REPLY == 0 )); then
    exit 1
  else
    if (( REPLY <= -1 || REPLY >= (1 + ${#options[@]}) )) ; then
      echo "$REPLY is invalid number..."
    else
      tocont=1
    fi
  fi
  echo
  if [[ ${tocont} == 1 ]] ; then
    break
  fi
 done
done
echo  "You picked number: $REPLY for ${opt}"
IDL_DM=$IDL_DLM
if [[ ${opt} =~ "linux.x86_64" ]]; then
    IDL_SO=$IDL_Linux64
#elif [[ ${opt} =~ "linux.x86" ]]; then
#    IDL_SO=$IDL_Linux32
#    I32="-32"
# elif [[ ${opt} =~ sparc64 ]]; then
#   IDL_SO=$IDL_Sparc64
# elif [[ ${opt} =~ sparc ]]; then
#   IDL_SO=$IDL_Sparc32
#   I32="-32"
else
    echo 'Error: require IDL/ENVI bin path for Linux x86_64' >&2;
    exit 1;
fi
today=`date '+%Y%m%d'`
arr=(${opt//// })
for fn in ${arr[@]}; do
    if [[ ${fn} =~ idl.*[0-9]+ ]]; then
        selected=$fn
        break
    fi
done
TEMPDIR=$(mktemp -d ~/idl_save_${selected}_${today}.XXXXXXXX)
echo " "
echo "Backing up the old files and downloading the new ones to ${TEMPDIR}"
cd ${TEMPDIR}
echo " "
echo "$GET ${IDL_SO}"
$GET ${IDL_SO}
if [ $? -ne 0 ]; then
     echo "Error: failed $GET ${IDL_SO}" >&2
     exit 1
fi
echo " "
echo "$GET ${IDL_DM}"
$GET ${IDL_DM}
if [ $? -ne 0 ]; then
     echo "Error: failed $GET ${IDL_DM}" >&2
     exit 1
fi
echo " "
echo "$GET ${LEAPS}"
$GET ${LEAPS}
if [ $? -ne 0 ]; then
     echo "Error: failed $GET ${LEAPS}" >&2
     exit 1
fi
echo " "
echo "$GET ${SHOWVERS}"
$GET ${SHOWVERS}
if [ $? -ne 0 ]; then 
     echo "Error: failed $GET ${SHOWVERS}" >&2
     exit 1
fi
echo "$GET ${SHOWVERS2}"
$GET ${SHOWVERS2}
if [ $? -ne 0 ]; then 
     echo "Error: failed $GET ${SHOWVERS2}" >&2
     exit 1
fi
cd ${opt}
echo " "
date_time=`date +%Y-%m-%d`
echo "Enter administrator account password: "
if [ ! -f idl_cdf.so-orig-$date_time ]; then
    sudo mv -v idl_cdf.so idl_cdf.so-orig-$date_time
    if [ $? -ne 0 ]; then
      exit 1
    fi
    sudo mv -v idl_cdf.dlm idl_cdf.dlm-orig-$date_time
    if [ $? -ne 0 ]; then
      exit 1
    fi
fi
sudo mv -v ${TEMPDIR}/idl_cdf.so idl_cdf.so
if [ $? -ne 0 ]; then
    exit 1
fi
sudo mv -v ${TEMPDIR}/idl_cdf.dlm idl_cdf.dlm
if [ $? -ne 0 ]; then
    exit 1
fi
cd ${TEMPDIR}
echo " "
echo "Testing the installation... version: $VERSION"
echo " "
if [ -d ${opt}/../../resource/xprinter ]; then
    export IDL_DIR=${opt}/../.. && export LD_LIBRARY_PATH=.:${opt} && export XPPATH=${opt}/../../resource/xprinter && echo ".run showcdfversion2" | ${opt}/idl ${I32}
  else
    export IDL_DIR=${opt}/../.. && export LD_LIBRARY_PATH=.:${opt} && echo ".run showcdfversion2" | ${opt}/idl ${I32}
fi
echo " "
echo "The patch (for .so and .dlm) is installed in ${opt}"
echo " "
echo "NOTE: The original version (.dlm and .so) there is saved and renamed"
echo "      with \"-orig-$date_time\" appended."
echo " "
rm showcdfversion.pro
rm showcdfversion2.pro

