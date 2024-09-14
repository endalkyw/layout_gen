class P:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def update(self, x_new, y_new):
        self.x = x_new
        self.y = y_new

    def transform_p(self, x_offset, y_offset):
        self.x = self.x + x_offset
        self.y = self.y + y_offset

    def to_tuple(self):
        t = (self.x, self.y)
        return t


