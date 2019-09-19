#!/usr/bin/env python3
'''
Commonly used functions for REST client applications.

* JSON and XML handled, parsed into objects.
* HTTP headers and POST data handled.

Jeremy Yang
10 Nov 2017
'''
import sys,os,re,time
import urllib,urllib.request,urllib.parse
import base64
import json
from xml.etree import ElementTree
from xml.parsers import expat
#
REST_TIMEOUT=10
REST_RETRY_NMAX=10
REST_RETRY_WAIT=5
#
##############################################################################
def GetURL(url,headers={},parse_json=False,usr=None,pw=None,parse_xml=False,nmax_retry=REST_RETRY_NMAX,verbose=0):
  '''Entry point for GET requests.'''
  return RequestURL(url,headers,None,usr,pw,parse_json,parse_xml,nmax_retry,verbose)

##############################################################################
def PostURL(url,headers={},data={},usr=None,pw=None,parse_json=False,parse_xml=False,nmax_retry=REST_RETRY_NMAX,verbose=0):
  '''Entry point for POST requests.'''
  return RequestURL(url,headers,data,usr,pw,parse_json,parse_xml,nmax_retry,verbose)

##############################################################################
def RequestURL(url,headers,data,usr,pw,parse_json,parse_xml,nmax_retry,verbose):
  '''Use internally, not externally.  Called by GetURL() and PostURL().  Only Basic authentication supported.'''
  if isinstance(data,dict):
    data=urllib.urlencode(data)
  req=urllib.request.Request(url=url,headers=headers,data=data)
  if usr and pw:
    req.add_header("Authorization", "Basic %s"%base64.encodestring('%s:%s'%(usr,pw)).replace('\n','')) 

  if verbose>1:
    print('url="%s"'%url,file=sys.stderr)
  if verbose>2:
    print('request type = %s'%req.type,file=sys.stderr)
    print('request method = %s'%req.get_method(),file=sys.stderr)
    print('request host = %s'%req.host,file=sys.stderr)
    print('request full_url = %s'%req.full_url,file=sys.stderr)
    print('request header_items = %s'%req.header_items(),file=sys.stderr)
    if data:
      print('request data = %s'%req.data,file=sys.stderr)
      print('request data size = %s'%len(req.data),file=sys.stderr)

  i_try=0
  while True:
    i_try+=1
    try:
      response=urllib.request.urlopen(req,timeout=REST_TIMEOUT)
      if not response:
        print('ERROR: urlopen failed.',file=sys.stderr)
        return None
      status=response.getcode()
      #print('DEBUG: status: {0}'.format(status), file=sys.stderr)
      #info=response.info()
      #url=response.geturl()
      fbytes=response.read()
      response.close()
    except urllib.request.HTTPError as e:
      if e.code==404:
        return ([])
      print('HTTP Error: %s'%e,file=sys.stderr)
      return None
    except urllib.request.URLError as e:  ## may be socket.error
      # may be "<urlopen error [Errno 110] Connection timed out>"
      print('URLError [try %d/%d]: %s'%(i_try,nmax_retry,e),file=sys.stderr)
      if i_try<nmax_retry:
        time.sleep(REST_RETRY_WAIT)
        continue
      return None
    except Exception as e:
      print('Error (Exception): %s'%e,file=sys.stderr)
      if i_try<nmax_retry:
        time.sleep(REST_RETRY_WAIT)
        continue
      return None
    break

  if fbytes.strip()=='': return None

  if parse_json:
    try:
      rval=json.loads(fbytes.decode('utf_8'), encoding='utf_8')
    except ValueError as e:
      try:
        ### Should not be necessary.  Backslash escapes allowed in JSON.
        ftxt_fix=fbytes.decode('utf_8').replace(r'\"','&quot;').replace(r'\\','')
        ftxt_fix=ftxt_fix.replace(r'\r','\\\\r') #ok?
        ftxt_fix=ftxt_fix.replace(r'\n','\\\\n') #ok?
        rval=json.loads(ftxt_fix, encoding='utf_8')
        if verbose>1:
          print('Apparently fixed JSON Error: %s'%e,file=sys.stderr)
        if verbose>2:
          print('DEBUG: Apparently fixed JSON: "%s"'%ftxt,file=sys.stderr)
          #sys.exit()
      except ValueError as e:
        if verbose:
          print('Failed to fix JSON Error: %s'%e,file=sys.stderr)
          print('DEBUG: ftxt_fix="%s"'%ftxt_fix,file=sys.stderr)
        rval=fbytes.decode('utf_8')
  elif parse_xml:
    try:
      rval=ParseXml(fbytes.decode('utf_8'))
    except Exception as e:
      if verbose: print('XML Error: %s'%e,file=sys.stderr)
      rval=fbytes.decode('utf_8')
  else: #raw
    rval=fbytes.decode('utf_8')

  return rval

#############################################################################
def ParseXml(xml_str):
  doc=None
  try:
    doc=ElementTree.fromstring(xml_str)
  except Exception as e:
    print('XML parse error: %s'%e,file=sys.stderr)
    return False
  return doc

#############################################################################
