import math
import sys

class Notes(object):
    def __init__(self):
        self.messages = []

    def addNotes(self, newMessage):
        self.messages.append(newMessage)
        return 0

    def printNotes(self):
        for i in self.messages :
            print(i)
        return 0

    def logNotes(self):
        notes = 'Start messages here\n'
        for i in self.messages :
            notes = notes + '\t' + i + '\n'
        notes = notes + 'End messages here\n'
        return notes
