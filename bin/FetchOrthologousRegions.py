
# coding: utf-8

# In[95]:

from __future__ import print_function

import urllib2
import BeautifulSoup as BS
from itertools import count
import re


# In[158]:

FB_BASE_URL="http://www.flybase.org"

def request_fb_gene(gene_name):
    query_template= FB_BASE_URL + "/cgi-bin/uniq.html?species=Dmel&cs=yes&db=fbgn&caller=genejump&context=%s"
    query = query_template % gene_name
    return query

def request_extended_gene_region(source, fbid, chromo):
    html_template = FB_BASE_URL + ""
    #source is species
    #id is flybase gene id
    #chromosome
    FASTA_template = FB_BASE_URL + "/cgi-bin/getseq.html?source=%s&id=%s&chr=%s&dump=PrecompiledFasta&targetset=gene_extended2000"
    the_url = FASTA_template % (source.lower(), fbid, chromo)
    #return the_url
    return urllib2.urlopen(the_url).read()


def request_chromo_region(name, chromo, start, stop):
	url_template = "http://flybase.org/cgi-bin/gbrowse2/%s/?plugin=FastaDumper;plugin_action=Go;name=%s:%i..%i"
	data = urllib2.urlopen(url_template % (name.lower(), chromo, start, stop) ).read()
	getter = re.compile("pre\\>(.*)\\</pre",flags=re.DOTALL)
	#import pdb
	#pdb.set_trace()
	foo = getter.findall(data)[0]
	return foo
# In[19]:




# In[20]:




# In[25]:




# In[42]:




# In[146]:

#FIND the ORTHOLOGS from the main search results.
def get_gene_orhologs(gene_name):
    foo = urllib2.urlopen(request_fb_gene(gene_name))
    text = foo.read()
    mydoc = BS.BeautifulSoup(text)
    children = mydoc.findChildren()

    def one_div_to_tuple(divdiv):
        chl = divdiv.findChildren("div")
        name = chl[0].text
        link = chl[2].findChildren("a")[0]
        return (name, link.text.split("\\")[0],
                link["href"].replace("/reports/","").replace(".html",""),
                link["href"], link.text)


    ortho_header = mydoc.find(id="gn15_1_1")
    all_ortho = ortho_header.nextSibling()[0]
    table_data = all_ortho.childGenerator()
    retlist = list()
    
    for one_row,row_num in zip(table_data,count()):
        #Skip first row
        if row_num < 1:
            first_pass = False
            continue

        retlist.append( one_div_to_tuple(one_row))
    return retlist


# In[135]:

#GET a gene's extended region when you know its URL
def get_gene_region(gene_ID):
    url = FB_BASE_URL + "/reports/" + gene_ID + ".html"
    data = urllib2.urlopen(url).read()
    data = data.replace("\n","")
    finder = re.compile("<th.*?Sequence location.*?</th><td>(.*)(?=</td>)")
    region = finder.findall(data)[0].split()[0]
    chromosome,other = region.split(":")
    start,stop = map(int,map(lambda x:x.replace(",",""),other.split("..")))
    return (chromosome, start, stop)


# In[173]:

def get_all_orthologs(gene_name,spacing=5000):
    ortho_names_ids = get_gene_orhologs(gene_name)
    #regions = [get_gene_region(i[2]) for i in ortho_names_ids]
    #print(ortho_names_ids)
    
    for the_name in ortho_names_ids:
        the_region = get_gene_region(the_name[2])
        #thefasta = request_extended_gene_region(the_name[1].lower(),the_name[2],the_region[0])
	direction = the_region[1] <= the_region[2]
	
	start = the_region[1] - spacing if direction else the_region[2] - spacing
	end = the_region[2] + spacing if direction else the_region[1] + spacing
	
        thefasta = request_chromo_region(the_name[1].lower(), the_region[0], start, end )
	yield (the_name[1].lower(), thefasta)
    


# In[172]:

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--space",type=int,default=5000,help="Number of nucleotides before and after the gene to take.")
parser.add_argument("GENE_NAME",type=str)
parser.add_argument("OUT_PREFIX",type=str)
args, other = parser.parse_known_args()

for one_name, fasta in get_all_orthologs(args.GENE_NAME, args.space ):
	with open(args.OUT_PREFIX+one_name+".fa","w") as outfile:
		outfile.write(fasta)


# In[ ]:



