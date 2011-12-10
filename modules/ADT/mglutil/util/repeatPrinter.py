import os

class RepeatPrinter:
    """
    The class allows printing repeating messages on the last line of a terminal (Unix only). This can be used to display a counter that does not scroll the terminal

    example:
        printer = RepeatPrinter()
        printer.prompt('Hello World')
        for i in range(20):
            printer.update(str(i))

    thsi code will display 'Hello World' followed by a counter
    """

    def __init__(self):
        # use tcup lines to find out how many lines in the terminal
        self.nbLines = int(os.popen("tput lines").readlines()[0])
        self.repeatingPosition = None


    def prompt(self, prompt):
        # print a message leading the repeating messages

        #os.system("tput sc")
        os.system("tput cup %d 0"%(self.nbLines-1))
        self.repeatingPosition = len(prompt) + 1  
        print prompt
        
    def update(self, msg):
        os.system("tput cup %d %d"%(self.nbLines-2, self.repeatingPosition))
        print msg



if __name__ == '__main__':
    from time import sleep

    printer = RepeatPrinter()
    printer.prompt('Hello World')

    sleep(0.5)
    for i in range(20):
        printer.update(str(i))
        sleep(0.2)

