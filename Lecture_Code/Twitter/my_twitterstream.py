
#
# Based on the excellent UW Coursera Data Science Course Code:
# https://github.com/uwescience/datasci_course_materials/blob/master/assignment1/twitterstream.py
#

import oauth2 as oauth
import urllib2 as urllib

###########################################################################
## Fill these in using the README instructions
access_token_key = "2236044578-kRh3NN4WFjngWwHal0a4WBCtSV4Mwk5phVinoEg"
access_token_secret = "NFV1BVb6WF09ZIAWHL1j5WFBGfahcvp48W4GMXb5tiFdu"
consumer_key = "bR288qpjcdJn3aNxsVzCA"
consumer_secret = "9I4dZrzMeOiCvJKBGmcscEn4DPDF2AJcHFMOVBhqumY"
###########################################################################

_debug = 0

oauth_token    = oauth.Token(key=access_token_key, secret=access_token_secret)
oauth_consumer = oauth.Consumer(key=consumer_key, secret=consumer_secret)

signature_method_hmac_sha1 = oauth.SignatureMethod_HMAC_SHA1()

http_method = "GET"


http_handler  = urllib.HTTPHandler(debuglevel=_debug)
https_handler = urllib.HTTPSHandler(debuglevel=_debug)

'''
Construct, sign, and open a twitter request
using the hard-coded credentials above.
'''
def twitterreq(url, method, parameters):
  req = oauth.Request.from_consumer_and_token(oauth_consumer,
                                             token=oauth_token,
                                             http_method=http_method,
                                             http_url=url,
                                             parameters=parameters)

  req.sign_request(signature_method_hmac_sha1, oauth_consumer, oauth_token)

  headers = req.to_header()

  if http_method == "POST":
    encoded_post_data = req.to_postdata()
  else:
    encoded_post_data = None
    url = req.to_url()

  opener = urllib.OpenerDirector()
  opener.add_handler(http_handler)
  opener.add_handler(https_handler)

  response = opener.open(url, encoded_post_data)

  return response

def fetchsamples():
  url = "https://stream.twitter.com/1/statuses/sample.json"
  parameters = []
  response = twitterreq(url, "GET", parameters)
  for line in response:
    print line.strip()

if __name__ == '__main__':
  fetchsamples()
