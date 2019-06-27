import pytest

from ornl.settings import unique_workspace_name as uwn
from ornl.sans.sns.eqsans.load import load_events


def test_load_events():
    ws_test_load_events= load_events('EQSANS_92353')
    assert ws_test_load_events.name() == 'ws_test_load_events'

    ws_name = uwn()
    ws = load_events('EQSANS_92353', output_workspace=ws_name)
    assert ws.name() == ws_name


if __name__ == '__main__':
    pytest.main()
