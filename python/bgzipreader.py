import struct as st
import argparse
import sys

def getoptions():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', required=True, action='store', dest='filename',
                            help='gzip file')
    return arg_parser.parse_args()


def main():
    args =getoptions()
    #getBamGzipBlockStats(args.filename)
    #splitblocks(args.filename)
    outputemptyblock(args.filename)
    
class bamfeild:
    types = {1:'B',2:'H',4:'I'}
    sizes = {'ID1':1,'ID2':1, 'CM':1, 'FLG':1, 'MTIME':4, 'XFL':1,'OS':1, 
             'XLEN':2,'SI1':1,'SI2':1,'SLEN':2,'BSIZE':2, 'CDATA':1, 'CRC32':4, 'ISIZE':4} 
    def __init__(self,bytes,name,number):
        self.bytes = bytes
        self.name= name
        self.number = number


    def unpack(self):
        return st.unpack("<{}{}".format(self.number,self.types[self.sizes[self.name]])
                         ,self.bytes)[0]


class bamfilecompressed:

    
    order = ['ID1','ID2', 'CM', 'FLG', 'MTIME', 'XFL','OS', 'XLEN','SI1','SI2', 'SLEN',
            'BSIZE', 'CDATA', 'CRC32', 'ISIZE'] 

    def __init__(self, ID1,ID2, CM, FLG, MTIME, XFL,OS, XLEN,SI1,SI2,
                 SLEN,BSIZE, CDATA, CRC32, ISIZE):
        self.ID1 = ID1
        self.ID2 = ID2
        self.CM = CM
        self.FLG =FLG
        self.MTIME = MTIME
        self.XFL = XFL
        self.OS = OS
        self.XLEN = XLEN
        self.SI1 = SI1
        self.SI2 = SI2
        self.SLEN = SLEN
        self.BSIZE = BSIZE 
        self.CDATA = CDATA
        self.CRC32 = CRC32
        self.ISIZE = ISIZE
        self.is_emp = self.is_empty()

    @classmethod
    def empty(cls):
        ID1 = bamfeild(st.pack('B',31),'ID1',1)
        ID2 = bamfeild(st.pack('B',139),'ID2',1)
        CM = bamfeild(st.pack('B',8), 'CM', 1)
        FLG = bamfeild(st.pack('B',4),'FLG',1)
        MTIME = bamfeild(st.pack('<I',0),'MTIME',1)
        XFL = bamfeild(st.pack('B',0),'XFL',1)
        OS = bamfeild(st.pack('B',255),'OS',1)
        XLEN = bamfeild(st.pack('<H',6),'XLEN',1)
        SI1 = bamfeild(st.pack('B',66),'SI1',1)
        SI2 = bamfeild(st.pack('B',67),'SI2',1)
        SLEN = bamfeild(st.pack('<H',2),'SLEN',1)
        BSIZE = bamfeild(st.pack('<H',27),'BSIZE',1)
        CDATA = bamfeild(st.pack('2B',3,0),'CDATA',1)
        CRC32 = bamfeild(st.pack('<I',0),'CRC32',1)
        ISIZE = bamfeild(st.pack('<I',0),'ISIZE',1)
        return cls(ID1,ID2, CM, FLG, MTIME, XFL,OS, XLEN,SI1,SI2, SLEN,
                            BSIZE, CDATA, CRC32, ISIZE)
        
    @classmethod
    def fromfile(cls, filehandle):
        ID1 = bamfeild(filehandle.read(1),'ID1',1)
        ID2 = bamfeild(filehandle.read(1),'ID2',1)
        CM = bamfeild(filehandle.read(1), 'CM', 1)
        FLG = bamfeild(filehandle.read(1),'FLG',1)
        MTIME = bamfeild(filehandle.read(4),'MTIME',1)
        XFL = bamfeild(filehandle.read(1),'XFL',1)
        OS = bamfeild(filehandle.read(1),'OS',1)
        XLEN = bamfeild(filehandle.read(2),'XLEN',1)
        SI1 = bamfeild(filehandle.read(1),'SI1',1)
        SI2 = bamfeild(filehandle.read(1),'SI2',1)
        SLEN = bamfeild(filehandle.read(2),'SLEN',1)
        BSIZE = bamfeild(filehandle.read(2),'BSIZE',1)
        cdatalen= BSIZE.unpack()-XLEN.unpack()-19
        CDATA = bamfeild(filehandle.read(cdatalen),'CDATA',cdatalen)
        CRC32 = bamfeild(filehandle.read(4),'CRC32',1)
        ISIZE = bamfeild(filehandle.read(4),'ISIZE',1)
        return cls(ID1,ID2, CM, FLG, MTIME, XFL,OS, XLEN,SI1,SI2, SLEN,
                    BSIZE, CDATA, CRC32, ISIZE)

    def is_empty(self):
        return self.bytes() == b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
        
    def get_block_stats(self):
        return "{}, {}, {}\n".format(self.ISIZE.unpack(), self.CRC32.unpack(), self.BSIZE.unpack())


    def items(self):
        return (self.ID1,self.ID2,self.CM,self.FLG,self.MTIME,self.XFL,self.OS,self.XLEN,
               self.SI1,self.SI2,self.SLEN,self.BSIZE,self.CDATA,self.CRC32,self.ISIZE)

    def bytes(self):
        b = b''
        for i in self.items():
            b += i.bytes
        return b
        
    def unpacked(self):
        tup = []
        for i in self.items():
            tup.append([i.name, i.unpack()])
        return tup
        

def splitblocks(filename):
    emptyBlock = bamfilecompressed.empty()
    with open(filename) as fh:
        files = 0
        startnew = True
        while(True):
            if startnew:
                files += 1
                outh = open("{:0=4d}.bam".format(files),'wb')
                startnew = False
            a = bamfilecompressed.fromfile(fh)
            outh.write(a.bytes())
            if a.ISIZE.unpack() != 65494:
                startnew = True
                outh.close()                
            if a.is_empty():
                outh.close()
                break


def getBamGzipBlockStats(filename):
    with open(filename) as fh:
        with open("blockstats.csv",'w') as out:
            while(True):
                block = bamfilecompressed.fromfile(fh)
                out.write(block.get_block_stats())
                if block.is_empty():
                    print(block.unpacked())
                    print(block.bytes())
                    break
                if block.ISIZE.unpack()==0:
                    print(block.unpacked())
                    break
            
            
def outputemptyblock(filename):
    emptyBlock = bamfilecompressed.empty()
    with open(filename, 'wb') as fh:
        fh.write(emptyBlock.bytes())




def nextfilehandle(nfile,stat,outh):
    pass

def read_gzip_record(fh,number,outh):
    pass
    
    
    
    
        
        

if __name__=='__main__':
    main()
