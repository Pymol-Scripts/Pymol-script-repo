from Example_services import EchoResponseWrapper

def echo(message):
    response = EchoResponseWrapper()
    response._Message = message
    return response
