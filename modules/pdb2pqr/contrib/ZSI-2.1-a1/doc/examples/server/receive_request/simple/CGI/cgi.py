from Echo_services import EchoResponseWrapper

def echo(message):
    response = EchoResponseWrapper()
    response._Message = message
    return response


from ZSI import dispatch
dispatch.AsCGI()
