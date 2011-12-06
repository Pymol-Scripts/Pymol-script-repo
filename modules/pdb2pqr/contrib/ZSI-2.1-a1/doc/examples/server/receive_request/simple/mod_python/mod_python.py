from ZSI import dispatch
from mod_python import apache

mod = __import__('encodings.utf_8', globals(), locals(), '*')
mod = __import__('encodings.utf_16_be', globals(), locals(), '*')

import MyHandler
def handler(req):
    dispatch.AsHandler(modules=(MyHandler,), request=req)
    return apache.OK
