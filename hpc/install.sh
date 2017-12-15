# $UWHPSC/notes/install.sh
# 
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install xfce4
sudo apt-get install jockey-gtk
sudo apt-get install xdm
sudo apt-get install ipython
sudo apt-get install python-numpy
sudo apt-get install python-scipy
sudo apt-get install python-matplotlib
sudo apt-get install python-dev
sudo apt-get install git
sudo apt-get install python-sphinx
sudo apt-get install gfortran
sudo apt-get install openmpi-bin
sudo apt-get install liblapack-dev
sudo apt-get install thunar
sudo apt-get install xfce4-terminal

# some packages not installed on the VM 
# that you might want to add:

sudo apt-get install gitk               # to view git history
sudo apt-get install xxdiff             # to compare two files
sudo apt-get install python-sympy       # symbolic python
sudo apt-get install imagemagick        # so you can "display plot.png"


sudo apt-get install python-setuptools  # so easy_install is available
sudo easy_install nose                  # unit testing framework
sudo easy_install StarCluster           # to help manage clusters on AWS
