import ConfigParser
config = ConfigParser.ConfigParser()
config.read("settings.ini")

config.get("Templates", "NAMES_TEMPLATE")