from sqlobject import *
import os.path, sys
from Bio import SeqIO
from Bio.Blast import NCBIXML

dbfile = 'blast_1.db3'

def init(new=False):
    conn_str = os.path.abspath(dbfile)
    conn_str = 'sqlite:'+ conn_str
    sqlhub.processConnection = connectionForURI(conn_str)
    if new:
        Proteome.dropTable(ifExists=True)
        Align.dropTable(ifExists=True)
        Proteome.createTable()
        Align.createTable()

class Proteome(SQLObject):
    pid = IntCol(alternateID=True)
    pdesc = StringCol()
    pseq = StringCol()
    aligns = MultipleJoin("Align")

class Align(SQLObject):
    proteome=ForeignKey("Proteome")
    querydesc = StringCol()
    subid=IntCol()
    subjectdesc= StringCol()
    evall = FloatCol()
    score=IntCol()
    queryseq=StringCol()
    match=StringCol()
    subjctseq= StringCol()
    align_length=IntCol()
    Identity=IntCol()
    ident_percent=FloatCol()
    
    
def load_proteome_table(file1,file2):
    files=[file1,file2]
    descs = []
    seqs = []
    for i in files:
        lines=(open(i, 'r').read().split('\n'))   
        data = ''
        for line in lines:
            if line.startswith('>'):
                if data:   
                    seqs.append(data)
                    data = ''
                descs.append(line)
            else:
                data += line.rstrip('\r\n')
        seqs.append(data)
    init(new=True)
    for i,j in zip(descs,seqs):
        a=i.split('|')
        t = Proteome(pid=int(a[1]), pdesc=i,pseq=j)

def load_align_table(file1,file2):
    init()
    files=[file1,file2]
    for i in files:
        result_handle = open(i)
        for blast_result in NCBIXML.parse(result_handle):
            for alignment in blast_result.alignments:
                for hsp in alignment.hsps:
                    title=alignment.title
                    e_val=hsp.expect
                    scr=hsp.score
                    a_len=alignment.length
                    ident=hsp.identities
                    q_seq=hsp.query
                    m_seq=hsp.match
                    s_seq=hsp.sbjct
                    q_row=blast_result.query.split('|')
                    pid=int(q_row[1])
                    r=Proteome.byPid(pid)
                    s_row=title.split('|')
                    
                    t=Align(querydesc=blast_result.query,
                            subid=int(s_row[3]),
                            subjectdesc=title,
                            evall=e_val,
                            score=scr,
                            align_length=a_len,
                            Identity=ident,
                            queryseq=q_seq,
                            match=m_seq,
                            subjctseq=s_seq,
                            ident_percent=(ident/float(a_len)),
                            proteome=r)

def protein_id():
    init()
    hs_name= Proteome.selectBy()
    for i in hs_name:
        print 'protein id: ', i.pid
        print 'protein desc: ',i.pdesc
        

def blast_result(p_id):
    init()

    try:
        hs_name = Proteome.selectBy(pid=int(p_id))[0]
        t= hs_name.aligns[0]
        print "query: ",t.querydesc
        print "subject id: ",t.subid
        print "subject: ",t.subjectdesc
        print "\nevalue: ",t.evall,
        print "\tscore: ",t.score
        print "alignment length: ",t.align_length
        print "identity: ",t.Identity
        print "identity percentage: ",t.ident_percent
        print "\n***Alignment***"
        qs= t.queryseq
        m= t.match
        ss= t.subjctseq
        print_align(qs,m,ss)
        return t.subid,t.align_length,t.Identity,t.evall
    except:
        print >>sys.stderr, "ID ",p_id,"does not exist"

    

def blast(p_id):
    a=blast_result(p_id)
    
def ortholog(p_id):
    a=blast_result(p_id)
    identity1=(a[2]/float(a[1]))
    if a is not None:
        print '\n\n'
        b=blast_result(a[0])
        identity2=(b[2]/float(b[1]))
        print '\n\n'
        if (b[0]==p_id) and identity1 >.6 and identity2 >.6 and a[3] < 1e-10 and b[3] < 1e-10:
            print "both proteins are Functional ORTHOLOG"
        else:
            print "both proteins are not Functional ORTHOLOG"
    
    
def print_align(qs,m,ss):
    term_width = 160
    all_lists = (qs, m, ss)
    length = max(map(len, all_lists))
    for offset in xrange(0, length, term_width):
        print
        for l in all_lists:
            print ''.join(''.join(map(str, l[offset:offset+term_width])))


        

