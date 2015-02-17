import os,subprocess,sys,glob,shutil

def getDirs(foldername):
    return os.walk(foldername).next()[1]

def Build(directory):
    cwd = os.getcwd()
    os.chdir(directory)
    p = subprocess.Popen('python setup.py build',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
    p.wait()
    if(p.returncode != 0 and p.returncode != None):
             print "     ----------------------------------------------------------"
             print "     > ERROR: MC-SPAM transit code couldn't be installed."
             print "     > Problem building code in "+directory+". The error was:\n"
             out, err = p.communicate()
             print spaced(err,"\t \t")
             print "     > If you can't solve the problem, please communicate"
             print "     > with the team for help.\n \n"
             os.chdir(cwd)
             sys.exit()
    libfolder = getDirs('build/.')
    for name in libfolder:
               if(name[0:3]=='lib'):
                 filename = glob.glob('build/'+name+'/*')
                 shutil.copy2(filename[0],'../../.')
    shutil.rmtree('build')
    os.chdir(cwd)

version = 'v.1.0'

print '\n'
print '\t ###############################################################'
print ''
print '\t                   MC-SPAM '+version+'\n'
print ''
print '\t         Authors: Nestor Espinoza (nespino@astro.puc.cl)'
print '\t                  Andres Jordan   (ajordan@astro.puc.cl)'
print ''
print '\t DISCLAIMER: If you use the results or parts of the code, please'
print '\t             cite Espinoza & Jordan (2015).'
print ''
print '\t > Installing transit codes...'
Build('main_codes/transit_code/planet/Python')
print '\t > ...done! Now you can safely use the code. See the README file \n'

