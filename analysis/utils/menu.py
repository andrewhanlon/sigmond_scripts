class Menu:

  def __init__(self, question, title=None):
    self.title = title
    self.question = question
    self.items = list()

  def addItem(self, item):
    self.items.append(item)

  def addItems(self, *items):
    self.items.extend(items)

  def getInput(self, multi=True, no_items=False):

    print()
    if self.title:
      print(self.title)

    item_menu = list()
    if multi:
      item_menu.append("[0] All")

    for item_count, item in enumerate(self.items, 1):
      item_menu.append(f"[{item_count}] {item}")

    if no_items:
      item_menu.append(f"[{item_count+1}] None")

    print("\n".join(item_menu))

    if multi:
      print("Specify choices separated by a space and/or ranges with a '-'")

    if len(self.items) == 1 and not no_items:
      print("{} {}".format(self.question, self.items[0]))
      return self.items if multi else self.items[0]
    else:
      answer = input("{} ".format(self.question)).split()

    if multi and "0" in answer:
      return self.items

    if no_items and "{item_count+1}" in answer:
      return None

    if not multi:
      if len(answer) != 1:
        print("Must provide one and only one answer")
        return self.getInput(multi=multi, no_items=no_items)

      try:
        return self.items[int(answer[0])-1]
      except ValueError:
        print("Invalid entries")
        return self.getInput(multi=multi, no_items=no_items)

    try:
      choices = list()
      for choice in answer:
        if '-' in choice:
          range_tup = choice.split('-')
          choices.extend(list(range(int(range_tup[0])-1, int(range_tup[1]))))
        else:
          choices.append(int(choice)-1)
      return [self.items[x] for x in choices]
    except (ValueError, IndexError):
      print("Invalid entries")
      return self.getInput(multi=multi, no_items=no_items)
