# -*- coding: utf-8 -*-
import sys
import urllib
import json
import shutil
import tempfile
from rdkit import Chem
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from concurrent import futures

def get_sd_file(aids, sdname,remove_inconclusive):
    with ThreadPoolExecutor(len(aids)) as executor:
        mols = set(executor.submit(get_sd_for_aid,aid,remove_inconclusive) for aid in aids)
    seen = set()
    count = 0
    writer = Chem.SDWriter(sdname)
    for future in futures.as_completed(mols):
        for m in future.result():
            if not m.HasProp("PUBCHEM_COMPOUND_CID"):
              continue
            cid =  m.GetProp("PUBCHEM_COMPOUND_CID")
            if cid in seen:
                continue
            seen.add(cid)
            writer.write(m)
            count += 1
    return count

def get_sd_for_aid(aid,remove_inconclusive):
        result = []
        cid_to_act = dict()
        aid = aid.strip()
        if not aid:
            return result
        csvfh = urllib.urlopen("http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+aid+"/CSV")
        header_line = csvfh.readline()
        if header_line.startswith("PUBCHEM_SID,PUBCHEM_CID,PUBCHEM_ACTIVITY_OUTCOME"):
            header_line = header_line.rstrip()
            header = header_line.split(",")
            for line in csvfh:
                line = line.rstrip();
                data = line.split(",")
                if len(data) != len(header):
                    continue
                cid_to_act[data[1]] = dict()
                for i,d in enumerate(data):
                    if i > 1:
                        cid_to_act[data[1]][header[i]] = d
        csvfh.close()
        listkeyfh = urllib.urlopen("http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+aid+"/cids/JSON?list_return=listkey")
        jsonlistkey = listkeyfh.read()
        listkeyfh.close()
        listkey = json.loads(jsonlistkey)
        if not listkey.has_key("IdentifierList"):
            return result
        if not listkey["IdentifierList"].has_key("ListKey"):
            return result
        key = listkey["IdentifierList"]["ListKey"]
        if not key:
            return result
        sdffh = urllib.urlopen("http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/"+key+"/SDF")
        if sdffh.getcode() != 200:
            return result
        tmpsdf = tempfile.NamedTemporaryFile()
        shutil.copyfileobj(sdffh,tmpsdf)
        sdffh.close()
        for m in Chem.SDMolSupplier(tmpsdf.name):
            if m is not None:
                if not m.HasProp("PUBCHEM_COMPOUND_CID"):
                  continue
                cid =  m.GetProp("PUBCHEM_COMPOUND_CID")
                if remove_inconclusive and cid_to_act[cid]["PUBCHEM_ACTIVITY_OUTCOME"] != "Active" and cid_to_act[cid]["PUBCHEM_ACTIVITY_OUTCOME"] != "Inactive":
                    continue
                m.SetProp("PUBCHEM_AID",str(aid))
                for k,v in cid_to_act[cid].items():
                    m.SetProp(str(k),str(v))
                result.append(m)            
        tmpsdf.close()
        return result


def get_aids(term,start=0):
    eufh = urllib.urlopen("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pcassay&term="+term+"&RetStart="+str(start))
    xmlstr = eufh.read()
    root = ET.fromstring(xmlstr)
    aids = [aid.text for aid in root.findall("./IdList/Id")]
    step = root.find("./RetMax").text;
    total = root.find("./Count").text;
    return step,total,aids

def get_description(aid):
 try:
  descfh = urllib.urlopen("http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+aid+"/summary/JSON")
  jsonstr = descfh.read();
  jsonobj = json.loads(jsonstr)
  name = jsonobj["AssaySummaries"]["AssaySummary"][0]["Name"]
  num = jsonobj["AssaySummaries"]["AssaySummary"][0]["CIDCountAll"]
 except KeyError:
  return None
 return (name,num,aid)

def get_list_aids(query, start=0):
    if query is None:
        choices = [] 
        return choices
    (step,total,aids) = get_aids(query,start)
    choices = []
    with ThreadPoolExecutor(len(aids)) as executor:
        future_choices = set(executor.submit(get_description,aid) for aid in aids)
    for future in futures.as_completed(future_choices):
     if future.result() is not None:
        (name,num,aid) = future.result()
        name += " ("+str(num)+")"
        choices.append((aid,name))
    choices.sort()
    return (step,total,choices)   

def get_url_base():
    return "http://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid="

def print_aid(aid):
    print aid, get_description(aid)

if __name__ == '__main__':
    query = sys.argv[1]
    sdname = sys.argv[2]
    remove_inconclusive = True

    (step,total,choices) = get_list_aids(query)
    # 20 2
    aids = [c[0] for c in choices]
    num_mols = get_sd_file(aids,sdname,remove_inconclusive)
    # 52s
    # 6s
