def Debugger(*args):
    for arg in args:
        try :
            print(float(arg))
            continue
        except ValueError:
            print(arg)
            continue
        except TypeError:
            print(arg)
            continue
        else :
            print()

    quit()
