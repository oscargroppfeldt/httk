#
#    The high-throughput toolkit (httk)
#    Copyright (C) 2012-2018 Rickard Armiento
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Do import inside class __init__ so that the missing import is only triggered if the class is actually used.
import os
import six

class TemplateEngineTemplator(object):
    def __init__(self, template_dir, template_filename, base_template_filename = None):
        try:
            from web.template import render
        except ImportError:
            raise Exception("Missing web.py module.")
        self.render = render

        self.template_dir = template_dir
        self.template_filename = template_filename
        self.template_name = template_filename[:-len(".html")] if template_filename.endswith(".html") else template_filename
        self.filename = os.path.join(template_dir, template_filename)

        self.dependency_filenames = [self.filename]
        if base_template_filename is not None:
            self.base_filename = os.path.join(self.template_dir,base_template_filename)
            self.dependency_filenames += [self.base_filename]
            self.base_template = base_template_filename.split(os.extsep)[0]
        else:
            self.base_filename = None

    def apply(self, content = None, data = None, *subcontent):
        if data is None:
            data = {}
        else:
            self.data = dict(data)
        if self.base_filename is not None:
            templator = self.render(self.template_dir,base=self.base_template,globals=data)
            output = six.text_type(getattr(templator,self.template_name)(content,*subcontent))
        else:
            templator = self.render(self.template_dir,globals=data)
            output = six.text_type(getattr(templator,self.template_name)(content,*subcontent))

        return output

    def get_dependency_filenames(self):
        return self.dependency_filenames
