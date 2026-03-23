import os
import lat_binned_analysis

install_dir = os.path.split(lat_binned_analysis.__file__)[0]

def main():

    # Copy starting files to new analysis directory:
    new_dir = os.getcwd()
    os.system("scp %s/inputs.yaml %s/client.py %s" %(install_dir,install_dir,new_dir))

########################
if __name__=="__main__":
        main(sys.argv)
