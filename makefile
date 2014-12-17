plugin: eclipse_reader_template.xml eclipse_reader_info.py eclipse_reader_data.py make_plugin_xml.py
	python make_plugin_xml.py eclipse_reader_template.xml eclipse_reader_info.py eclipse_reader_data.py
	
package: plugin
	rm -rf package
	mkdir package
	cp eclipse_reader.xml README.md Zprava_projektu.txt show_wells.py package
	cd package; \
	zip ../package *