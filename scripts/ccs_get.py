#!/usr/bin/env python
""" Get raw data for pacbio from Mount Sinai """

import sys , io ,os , urllib
import optparse, logging


class ReadGetter:
    def __init__( self ):
        self.__parseArgs( )
        self.__initLog( )


    def __parseArgs( self ):
        """Handle command line argument parsing"""
        
        usage = "%prog [--help] [options] JOBID"
        parser = optparse.OptionParser( usage=usage, description=__doc__ )

        parser.add_option( "-l", "--logFile", help="Specify a file to log to. Defaults to stderr." )
        parser.add_option( "-d", "--debug", action="store_true", help="Increases verbosity of logging" )
        parser.add_option( "-i", "--info", action="store_true", help="Display informative log entries" )
        parser.add_option( "-p", "--profile", action="store_true", help="Profile this script, dumping to <scriptname>.profile" )
        parser.add_option( "-t", "--prefix", help="prefix to append at beginning of the file")
        parser.add_option( "-e", "--ext", help="bas.h5, ccs.fasta, etc.")
        parser.add_option( "-s", "--server", help="Change server (default is http://node1.1425mad.mssm.edu/)")
        parser.add_option( "--noprefix", action="store_true", help="Don't put a prefix in front of the file name - just use whatever the name on the server is")
        parser.set_defaults( logFile=None, debug=False, info=False, profile=False, ext="bas.h5", prefix=False, server='http://node1.1425mad.mssm.edu/', noprefix = False)
        
        self.opts, args = parser.parse_args( )

        if len(args) != 1:
            parser.error( "Expected 1 argument." )
        self.server = self.opts.server    
        self.inFN = args[0]

        self.jobNumber = str(args[0])

        if self.opts.prefix == False:
            self.prefix = 'jobID_'+self.jobNumber+'_' 
        else:
            self.prefix = self.opts.prefix
        
    def __initLog( self ):
        """Sets up logging based on command line arguments. Allows for three levels of logging:
        logging.error( ): always emitted
        logging.info( ): emitted with --info or --debug
        logging.debug( ): only with --debug"""

        logLevel = logging.DEBUG if self.opts.debug else logging.INFO if self.opts.info else logging.ERROR
        logFormat = "%(asctime)s [%(levelname)s] %(message)s"
        if self.opts.logFile != None:
            logging.basicConfig( filename=self.opts.logFile, level=logLevel, format=logFormat )
        else:
            logging.basicConfig( stream=sys.stderr, level=logLevel, format=logFormat )
                                                                 
    def run( self ):
        """Executes the body of the script."""
    
        logging.info("Log level set to INFO")
        logging.debug("Log Level set to DEBUG")
        
        jobNumber = self.jobNumber

        #print jobNumber
        #print jobNumber[:3]

        portalURL = self.server
        home = '/home/sbsuser/pacbio/raw/'
        splicehome = '/pacbio/raw'
        ext = self.opts.ext
        #print ext
        records = set()
            
        if ext == "ccs_reads.fastq":
            cmd = 'wget http://node1.1425mad.mssm.edu/pacbio/secondary/%s/%s/data/ccs_reads.fastq' % (jobNumber[:3],jobNumber)
            logging.info(cmd)
            os.system(cmd)
        elif ext == "cmp.h5":
            #logIn = urllib.urlopen('http://node1.1425mad.mssm.edu/pacbio/secondary/%s/%s/data/aligned_reads.cmp.h5' % (jobNumber[:3],jobNumber))
            cmd = 'wget http://node1.1425mad.mssm.edu/pacbio/secondary/%s/%s/data/aligned_reads.cmp.h5' % (jobNumber[:3],jobNumber)
            logging.info(cmd)
            os.system(cmd)
        elif ext == "bax.h5":

            logIn = urllib.urlopen('http://node1.1425mad.mssm.edu/pacbio/secondary/%s/%s/log/smrtpipe.log' % (jobNumber[:3],jobNumber))
            for line in logIn:
                if home in line:
                    ll = line.split(" ")
                    #print "starting a new set:"
                    #print ll
                    #print "raw: ", line
                    line = ll[-3]
                    #print "split: ", line
                    #print portalURL
                    #print line[line.find(splicehome):line.find("bax.h5")]
                    #print ext
                    if not "m" in line[line.find(splicehome):line.find("bax.h5")]:
                        continue

                    records.add(portalURL+line[line.find(splicehome):line.find("bax.h5")] + ext)
                    records.add(portalURL+line[line.find(splicehome):line.find("bax.h5")-2] + "bas.h5")
                    #print records
        else:
            print >>sys.stderr, "Not supported file type!"

        for address in records:
            logging.info(address)
            fileIn = urllib.urlopen(address)
            if self.opts.noprefix:
                fileout = open(address.split('/')[-1],'w')
            else:
                fileout = open(self.prefix+address.split('/')[-1],'w')
            fileout.write(fileIn.read())
            fileout.close()	
        return 0

if __name__ == "__main__":
    app = ReadGetter()
    if app.opts.profile:
        import cProfile
        cProfile.run( 'app.run()', '%s.profile' % sys.argv[0] )
    sys.exit( app.run() )

#for address in records:
#    print address
#    filename = address.split('/')[-1]
#    os.system('wget %s >> %sjobID_%s_%s' % (address , savepath, jobNumber , filename))
