import sys, os

this_file = __file__
abs_path_root = '/'.join(this_file.split('/')[0:-1])
sys.path.append(abs_path_root)
os.chdir(abs_path_root)

import jinja2
import time 
import marshal
import re
import zipfile
import glob
import web
import pickle

from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline

from abiparse3 import abiparse
from UtilityFunctions import is_fasta

##

urls = (
      '/', 'FrontHandler',
      '/?','FrontHandler',
      '/abiupload','ABIUploadHandler',
      '/align','AlignHandler',
      '/upload_seq','SeqUploadHandler',
      '/upload','UniqueUploadHandler',
      '/trouble','TroubleHandler',
      '/export','ExportHandler',
      '/abiprocess','ABIProcessingHandler',
      '/process/?','ProcessingHandler',
      '/seqprocess','SeqProcessingHandler',
      '/serve/(.+)','ServeHandler',
      )




# silence database logging
# web.config.debug = False

class report:
    def __init__(self):
        self.num = None
        self.groups = None

    def setNum(self,num):
        self.num = num

    def setGroups(self,groups):
        self.groups = groups

class UniqueSeq():

    def __init__(self):
        self.db = web.database(dbn='sqlite',db=':memory:')
        self.db.query('''CREATE TABLE seq_list (
             seq_id TEXT,
             seq TEXT,
             seq_pickle TEXT,
             session TEXT);''')

        self.uniq_sequences = []
        self.seq_dict = {}
        self.start_ind = 0
        self.stop_ind = None
    
    def set_start_ind(self,value):
        self.start_ind = value
    def set_stop_ind(self,value):
        self.stop_ind = value
   
    def readLoop(self,parsed,session):
        for seq_record in parsed:
            #print seq_record.id
            seq_id_str = seq_record.id
            seq_str = seq_record.seq.tostring()
            seq_pickle = pickle.dumps(seq_record)
            session=session

            self.db.insert('seq_list',
                           seq_id=seq_id_str,
                           seq=seq_str[self.start_ind:self.stop_ind],
                           seq_pickle=seq_pickle,
                           session=session)

 
    def readFile(self,filename,filetype):
        session="imakey"
        
        handle = open(filename,'r') # will replace this with blob reader instance
        if filetype == 'fasta':
            parsed =  SeqIO.parse(handle,"fasta")
        elif filetype == 'genbank':
            parsed = SeqIO.parse(handle,"genbank")
        else: 
            # this check should not fail here
            # should already be checked
            raise web.seeother('/?err=type2')

        self.readLoop(parsed,session)


    def getUnique(self,session):


        # use a dictionary of lists, keyed to the unique sequence, containing the sequence ids
        # keep this dictionary as self.report, use write report to pickle and save this
        # can build more layers into this over time if other report elements become of interest

        entries = self.db.query("SELECT DISTINCT seq FROM seq_list WHERE session='%s';"%session)
        entries = list(entries)
        if not self.uniq_sequences:
            # compile list of all entries
            unique_entries = set([entry.seq for entry in entries])
            for entry in unique_entries:
                seq_query = self.db.query("SELECT * FROM seq_list WHERE seq='%s'"%entry)
                seq_query = list(seq_query)
                seq_obj = pickle.loads(str(seq_query[0].seq_pickle))
                self.uniq_sequences.append(seq_obj)

                # generate report of sequence coverage
                for seq_pick in seq_query:
                    seq_obj = pickle.loads(str(seq_pick.seq_pickle))
                    if self.seq_dict.has_key(seq_pick.seq):
                        self.seq_dict[seq_pick.seq].append(seq_obj.id)
                    else:
                        self.seq_dict[seq_pick.seq] = [seq_obj.id]
        

    def writeUnique(self,results_file,time_key):
        
        session = 'imakey' # will need to change this to time_key
        results_dir = '/uploads/results/'
 
        self.getUnique(session)

        output_handle = open(results_dir + time_key + results_file,'w')
        SeqIO.write(self.uniq_sequences, output_handle, "fasta")
        output_handle.close()

        return session 

    def writeReport(self,report_file,time_key):
        # want to reconfigure this so as to give not just the number of unique sequences, but
        # also more information about the set of sequences as a whole
        # start by being able to identify the sequence groups
        num = self.uniqueNumber('imakey')
        groups = self.getGroups('imakey')
        obj = report()
        obj.setNum(num)
        obj.setGroups(groups)

        file_dir = abs_path_root + '/uploads/reports/'

        pickle.dump(obj, open(file_dir + time_key + report_file, 'wb'))
    
    def getGroups(self,session):
        self.getUnique(session)
        return self.seq_dict
    
    def uniqueNumber(self,session):
        self.getUnique(session)
        return len(self.uniq_sequences)

## end seq picker stuff here

class Handler():
    def render_template(self,template_name, **context):
       extensions = context.pop('extensions', [])
       globals = context.pop('globals', {})

       jinja_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')),
            extensions=extensions,
            )
       jinja_env.globals.update(globals)

       #jinja_env.update_template_context(context)
       return jinja_env.get_template(template_name).render(context)

class FrontHandler(Handler):
    def GET(self):
        if not web.ctx.env['HTTP_USER_AGENT'].find('MSIE') == -1:
            return self.render_template('IE_error.html')

        query = web.input()
        if query:
           error = query.get('err',None)
           if error:
               if error == 'type':
                   error_message = "ERROR: The filetype you submitted is unsupported or you are using an unsupported browser. <br> <br> Currently only fasta and genbank formated files are accepted. <br><br> The site has been tested with Chrome 21.0 and Firefox 14.0. Current versions of these browsers should be compatible." 
               elif error == 'type2':
                   error_message = "ERROR: The filetype you submitted is unsupported.  This was caught after upload was validated and is hence a major error. <br><br> Please get in touch with Andrew"
               else:
                   error_message = "META-ERROR: The error submitted to handler is not a known error.<br><br> If you're suprsied to see this, please get in touch with Andrew."
           else:
               error_message = "Your query has non-recognized fields.<br><br> If you're surprised to see this, please get in touch with Andrew."
        else:
           error_message = ""
        return self.render_template('front_page.html',upload_url='/upload',error_message=error_message)
    
class TroubleHandler(Handler):
    def GET(self):
        return self.render_template('trouble_page.html')


class UniqueUploadHandler(Handler):
    def process_file(self,filename,time_key,start_ind,stop_ind,filetype):

        test = UniqueSeq()
        test.start_ind = int(start_ind)
        test.stop_ind = int(stop_ind)
        test.readFile(filename,filetype)

        # replace 'imakey' with time_key when you get a chance
        num = test.uniqueNumber('imakey')

        # write file here, return filename
        results_file = 'results.fasta'
        report_file  = 'report.p'

        test.writeUnique(results_file,time_key)
        test.writeReport(report_file,time_key)

        return time_key + results_file, time_key + report_file

    def POST(self):
        x = web.input(seqfile = {},start_ind=0,stop_ind=None)
        filedir = abs_path_root 
        if 'seqfile' in x:
            print x.seqfile.filename 
            time_key = str(time.time()).replace('.','')

            filepath=x.seqfile.filename.replace('\\','/')
            filename = time_key + filepath.split('/')[-1]
            
            filetype = filename.split('.')[-1]
            if filetype == 'fasta' or filetype == 'fas' or filetype == 'fna' or filetype == 'ffn' or filetype == 'faa' or filetype == 'frn':
                 filetype = 'fasta'
            elif filetype == 'gb':
                 filetype = 'genbank'
            else:
                 type_error = 'Unsuported file format'
                 print type_error, filepath
                 raise web.seeother('/?err=type')
            fout = open(filedir +'/uploads/'+ filename,'w')
            fout.write(x.seqfile.file.read())
            fout.close()

	results_file, report_file = self.process_file(filedir +'/uploads/'+ filename,
                                         time_key,
                                         x.start_ind,
                                         x.stop_ind,
                                         filetype)

        raise web.seeother('/process?results_file={0}&report_file={1}'.format(results_file,report_file))

class ABIUploadHandler(Handler):

    def POST(self):
        x = web.input(abifile = {}, cdrfile={})
        filedir = abs_path_root 
        #this if statement doesn't check for errors correctly
        if 'abifile' in x and 'cdrfile' in x:
            time_key = str(time.time()).replace('.','')
            filepath= x.abifile.filename.replace('\\','/')
            filename = time_key + 'abizipped.zip'

            filename_cdr = time_key + 'CDR.fas'


            print filepath
            #filetype = filename.split('.')[-1]
            #if filetype == 'zip':
            #    pass
            #else:
            #    raise web.seeother('/export?err=type')

            fout = open(filedir +'/uploads/'+ filename,'w')
            fout.write(x.abifile.file.read())
            fout.close()
            
            fout_cdr = open(filedir +'/uploads/'+ filename_cdr,'w')
            fout_cdr.write(x.cdrfile.file.read())
            fout_cdr.close()
            
            dirname = filedir + '/uploads/' + time_key + 'ABI'
            os.makedirs(dirname)

            abizip = zipfile.ZipFile(filedir +'/uploads/'+ filename,"r") 

            for info in abizip.infolist():
                 fname = info.filename
                 print fname
                 if fname.split('.')[-1] == 'ab1':
                     data = abizip.read(fname)
                     filename = dirname + '/'+ fname.split('/')[-1]
                     fout = open(filename,'w')
                     fout.write(data)
                     fout.close()
            abiparse(dirname + '/','*.ab1',filedir +'/uploads/'+ filename_cdr)

            zout = zipfile.ZipFile(filedir + '/uploads/results/' + time_key + 'abiexport.zip','w')
            for filename in glob.glob(dirname + '/*.*'):
                zout.write(filename,'abiexport/' + filename.split('/')[-1])
            zout.close() 

        raise web.seeother('/abiprocess?results_file={0}'.format(time_key + 'abiexport.zip'))


class ProcessingHandler(Handler):
    def getReport(self,report_file):
        file_dir = abs_path_root + '/uploads/reports/'
        print report_file
        report = pickle.load(open(file_dir +  report_file, 'rb'))
        return report.num, report.groups
    
    def GET(self):

        query = web.input()
        report_file = query.report_file
        file_name = query.results_file
        max_display_length = 70

        num, groups = self.getReport(report_file)
        if len(groups.keys()[0])>max_display_length:
            shorten = True
        else:
            shorten = False
        
        return self.render_template('processing_page.html',uniqueNumber=num,file_name=file_name,seq_groups=groups,shorten=shorten) 

class AlignHandler(Handler):
    def GET(self):
        if not web.ctx.env['HTTP_USER_AGENT'].find('MSIE') == -1:
            return self.render_template('IE_error.html')

        query = web.input()
        if query:
           error = query.get('err',None)
           if error:
               if error == 'type':
                   error_message = "ERROR: The filetype you submitted is unsupported or you are using an unsupported browser. <br> <br> Currently only fasta and genbank formated files are accepted. <br><br> The site has been tested with Chrome 21.0 and Firefox 14.0. Current versions of these browsers should be compatible."
               elif error == 'type2':
                   error_message = "ERROR: The filetype you submitted is unsupported.  This was caught after upload was validated and is hence a major error. <br><br> Please get in touch with Andrew"
               else:
                   error_message = "META-ERROR: The error submitted to handler is not a known error.<br><br> If you're suprsied to see this, please get in touch with Andrew."
           else:
               error_message = "Your query has non-recognized fields.<br><br> If you're surprised to see this, please get in touch with Andrew."
        else:
           error_message = ""
        return self.render_template('seq_page.html',upload_url='/upload_seq',error_message=error_message)
        

def align_file(filename,time_key,filetype):
    filedir = abs_path_root 
    results_file = 'align_results.fasta'

    cline = MuscleCommandline(input=filename, out=filedir + '/uploads/results/' +  time_key + results_file)
    cline()

    return time_key + results_file


class SeqUploadHandler(Handler):
    def align_file(self,filename,time_key,filetype):
        filedir = abs_path_root 
        results_file = 'align_results.fasta'

        cline = MuscleCommandline(input=filename, out=filedir + '/uploads/results/' +  time_key + results_file)
        cline()

        return time_key + results_file

    def POST(self):
        x = web.input(seqfile = {})
        filedir = abs_path_root 
        if 'seqfile' in x:
            time_key = str(time.time()).replace('.','')
            filepath= x.seqfile.filename.replace('\\','/')
            filename = time_key + filepath.split('/')[-1]
            
            print filepath
            filetype = filename.split('.')[-1]
            if filetype == 'zip':

                # handle the zipped file
                fout = open(filedir +'/uploads/'+ filename,'w')
                fout.write(x.seqfile.file.read())
                fout.close()

                dirname = filedir + '/uploads/' + time_key + 'FAS'
                os.makedirs(dirname)

                faszip = zipfile.ZipFile(filedir +'/uploads/'+ filename,"r")
                out_handle = open(filedir + '/uploads/results/' + time_key + 'multiseq.fasta','w')

                for info in faszip.infolist():
                     fname = info.filename
                     print fname
                     if is_fasta(fname): 
                         data = faszip.read(fname)
                         filename = dirname + '/'+ fname.split('/')[-1]
                         fout = open(filename,'w')
                         fout.write(data)
                         fout.close()
                         fin = open(filename,'r')
                         for record in SeqIO.parse(fin,"fasta"):
                             SeqIO.write(record,out_handle,"fasta")
                         fin.close()
                filename =  time_key + 'multiseq.fasta'
            else:
                time_key = str(time.time()).replace('.','')
    
                filepath=x.seqfile.filename.replace('\\','/')
                filename = time_key + filepath.split('/')[-1]

                filetype = filename.split('.')[-1]
                if filetype == 'fasta' or filetype == 'fas' or filetype == 'fna' or filetype == 'ffn' or filetype == 'faa' or filetype == 'frn':
                     filetype = 'fasta'
                else:
                     type_error = 'Unsuported file format'
                     print type_error, filepath
                     raise web.seeother('/align?err=type')
                fout = open(filedir +'/uploads/results/'+ filename,'w')
                fout.write(x.seqfile.file.read())
                fout.close()

        results_file = self.align_file(filedir +'/uploads/results/'+ filename,
                                         time_key,
                                         filetype)

        raise web.seeother('/seqprocess?results_file={0}'.format(results_file))

class ABIProcessingHandler(Handler):
    def GET(self):
        query = web.input()
        file_name = query.results_file
        return self.render_template('abi_processing_page.html',file_name=file_name)

class SeqProcessingHandler(Handler):
    def GET(self):
        query = web.input()
        file_name = query.results_file
        return self.render_template('seq_processing_page.html',file_name=file_name)

class ServeHandler():
    def GET(self,key):
        path = abs_path_root + '/uploads/results/' + key
        if os.path.exists( path ):
            getFile = file(path, 'r')
            web.header('Content-type','application/octet-stream')
            if path.split('.')[-1] == 'fasta':
                web.header("Content-Disposition","filename=results.fasta")
            elif path.split('.')[-1] == 'zip':
                web.header("Content-Disposition","filename=abiexport.zip")
            return getFile.read()
        else:
            raise web.notfound()

class ExportHandler(Handler):
    def GET(self):
        error=''
        return self.render_template('export_page.html',upload_url='/abiupload',error=error)

application = web.application(urls, globals()).wsgifunc()

