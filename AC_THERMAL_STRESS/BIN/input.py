def paramin(name):
    
    while True:
        try:
            x = float(input(name))
            break
        except ValueError:
            print('Oops!  That was no valid number.  Try again...')
    return x


def indexin(name):
   
    while True:
        try:
            x = int(input(name))
            break
        except ValueError:
            print('Please enter an integer!')
    return x