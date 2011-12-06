from ZSI import TC

class User:
    def __init__(self, name=None, userId=None, age=-1):
        self.Name = name
        self.UserId = userId
        self.Age = age

User.typecode = TC.Struct(User,
                          [TC.String('Name', oname='Name'),
                           TC.String('UserId', oname='UserId'),
                           TC.Iinteger('Age', oname='Age')],
                           'User')

