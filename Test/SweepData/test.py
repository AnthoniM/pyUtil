import numpy as np
import coreclasses as cc 

def js_iterator(nptsl):
    """ This increments the last element of js but
        carries to previous element any carryover above
        npts

        Taken from pyHegel.commands written by
        Christian Lupien
    """
    multiN = len(nptsl)
    js = [0]*multiN
    while True:
        yield js
        i = multiN -1
        while i >= 0:
            js[i] += 1
            if js[i] >= nptsl[i]:
                if i == 0:
                    return # we are done
                js[i] = 0
                i -= 1
            else:
                break

path="/home/toni/Codes/Python/modules/Test/SweepData/TestFiles/"

class Test(cc.SweepData):
    def __init__(self, field_dim, param_dim, ncols, nrows, head=True, base_filename="test", save_path="./TestFiles/", silent_fields=None, relevant_fields=None):
        self.field_dim=field_dim
        self.param_dim=param_dim
        self.ncols=ncols
        self.nrows=nrows
        self.head=head
        self.path=save_path
        self.relevant_fields=relevant_fields
        N = len("%s"%np.prod(self.param_dim))
        self.readval = "_readval_%0"+str(N)+"d.txt"
        self.labels = "abcdefghijklmnopqrstuvwxyz"
        self.base_filename = base_filename

        self.field_names=list(self.labels[:len(self.field_dim)])
        #self.pattern = self.base_filename+"_.*"+"_.*".join(["(%s\d+)"%i for i in self.field_names])+".*"
        self.comment = "#readback numpy shape for line part: "+", ".join([str(i) for i in self.param_dim])
        #if self.head:
        #    self.pattern+="readval_(\d+).*"
        #    self.field_names += ['_readval']
        #else:
        if not self.head:
            # Put the parameters in each file
            self.nrows=np.prod(self.param_dim)
            self.ncols+=len(self.param_dim)
        self.param_names=list(self.labels[len(self.field_dim):][:len(self.param_dim)])
        self.column_names=["col%s"%i for i in range(1,self.ncols+1)]
        print self.field_names, self.field_dim
        print self.param_names, self.param_dim
        print self.column_names
        self.generate_files()
        self.extra_dim=0
        if silent_fields:
            self.extra_dim=1
            print "Initial field_dim", self.field_dim, silent_fields
            mask = np.asanyarray(silent_fields).astype(bool)
            self.field_dim=list(np.asanyarray(field_dim)[mask])
            self.field_names=list(np.asanyarray(self.field_names)[mask])
        self.pattern = self.base_filename+"_.*"+"_.*".join(["(%s\d+)"%i for i in self.field_names])+".*"
        if self.head:
            self.pattern+="readval_(\d+).*"
            self.field_names += ['_readval']
        print self.pattern
        self.load()
        
    def generate_files(self):#self.field_dim, self.param_dim, self.ncols, self.nrows, self.head=True):


        base_filename = self.base_filename
        if self.head:
            for i in range(len(self.field_dim)):
                base_filename += "_%s"%self.labels[i]+"%s"
            base_filename += ".txt"
            for a in js_iterator(self.field_dim):
                self.head_data = ""
                self.head_filename = base_filename%tuple(a)
                self.head_data += '/'.join(["%s%s"%(i,j) for i,j in zip(self.labels[:len(a)],a)])
                with open(self.path+self.head_filename, 'w') as f:
                    f.writelines(self.comment+'\n')
                    for i,b in enumerate(js_iterator(self.param_dim)):
                        data = self.head_data
                        filename = self.head_filename[:-4]+self.readval%i
                        params = ["%s%s"%(j,k) for j,k in zip(self.labels[len(a):][:len(b)],b)]
                        data += '/'+'/'.join(params)
                        f.writelines('\t'.join(params)+'\n')
                        with open(self.path+filename, 'w') as g:
                            for r in range(self.nrows):
                                d = '\t'.join([data+"/col%s"%c+"/row%s/"%r for c in range(self.ncols)])
                                g.write(d+'\n')
        else:
            for i in range(len(self.field_dim)):
                base_filename += "_%s"%self.labels[i]+"%s"
            base_filename += ".txt"
            for a in js_iterator(self.field_dim):
                self.head_data = ""
                self.head_filename = base_filename%tuple(a)
                self.head_data += '/'.join(["%s%s"%(i,j) for i,j in zip(self.labels[:len(a)],a)])
                with open(self.path+self.head_filename, 'w') as g:
                    g.writelines(self.comment+'\n')
                    for i,b in enumerate(js_iterator(self.param_dim)):
                        data = self.head_data
                        params = ["%s%s"%(j,k) for j,k in zip(self.labels[len(a):][:len(b)],b)]
                        data += '/'+'/'.join(params)
                        d = '\t'.join([data+"/col%s"%c+"/row%s"%i for c in range(self.ncols)])
                        g.write(d+'\n')
    def load(self):
            ####################
            # Load Data
            ####################
            super(Test, self).__init__(self.pattern, self.path, param_names=self.param_names, field_names=self.field_names, verbose=True, relevant_fields=self.relevant_fields)

    def test(self):
        ####################
        # Test column_names
        ####################
        for i,c in enumerate(self.column_names):
            if not all(['col%s'%(i) in s for s in getattr(self, c).flatten()]):
                print getattr(self, c)
                raise Warning, "col%s does not pass the test !"%(i+1)
        print "Has the right column attributes !"

        ####################
        # Test rows
        ####################
        if self.nrows and self.head:
            tmp = np.rollaxis(self.data, -1)
            assert(len(tmp)==self.nrows)
            print "Good number of rows ! "

        ####################
        # Test field_names
        ####################
        tmp_field_names = self.field_names
        tmp_field_dim = self.field_dim
#        if self.relevant_fields:
#            mask = np.asanyarray(self.relevant_fields).astype(bool)
#            print mask
#            if len(mask)==len(self.field_names):
#                tmp_field_names = np.asanyarray(self.field_names)[mask]
#                tmp_field_dim = np.asanyarray(self.field_dim)[mask[:-1]]
#            else:
#                tmp_field_names = np.asanyarray(self.field_names)[:-1][mask]
#                tmp_field_dim = np.asanyarray(self.field_dim)[mask]
        field_names = [tmp_field_names[i] for i,d in enumerate(tmp_field_dim) if d>1]
        param_names = [self.param_names[i] for i,d in enumerate(self.param_dim) if d>1]
        print "Test for field_names with dimension higher then 1 : ", field_names
        print "Test for param_names with dimension higher then 1 : ", param_names
        for i,a in enumerate(field_names+param_names):
            b = np.rollaxis(self.data, i+1+int(self.extra_dim))
            for j,c in enumerate(b):
                if not np.all([a+str(j) in k for k in c.flatten()]):
                    raise Warning, "Bad representation of data !"
        print "Good representation of data !"

print "Test with readval"
test1 = Test([4,2,1],[5,6,1,2],2,10)
test1.test()

print "Test without readval"
test2 = Test([4,1,2],[5,6,1,1,1],2,10, head=False, base_filename="test2")
test2.test()

print "Test with readval and irrelevent fields"
test3 = Test([4,1,2,6,3,3],[5,6,1,1,1],2,10, head=True, base_filename="test3", silent_fields=[1,1,0,1,0,0])
test3.test()

print "Test without readval and irrelevent fields"
test4 = Test([4,1,2,6,2,3],[5,6,1,1,1],2,10, head=False, base_filename="test4", silent_fields=[1,1,0,1,0,0])
test4.test()
