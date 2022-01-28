import tables

class Template:
    def __init__(self, name):
        if name[-3:]!='.h5':
            name+='.h5'
        self.h5 = tables.open_file(name, 'r+')
        pass

    def create_group(self, name, shape):
        "For example sol000"
        return self

    def create_table(self, name, shape):
        "For example phase000"
        return self

    def create_array(self, name, shape):
        "For example val, weight, pol, ..."
        return self

if __name__ == '__main__':
    test = Template('test.h5')
    test.create_group


